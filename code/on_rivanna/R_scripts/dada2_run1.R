#install.packages("BiocManager", repos='http://cran.us.r-project.org') # WORKS
#BiocManager::install("ShortRead", dependencies = TRUE)

library("dada2")
library("ShortRead")
library("phyloseq")
library("ggplot2")
sink('/project/tumi/carey/crypto_micro/dada2_R_outfile_r1.txt')
packageVersion('dada2')

# get raw sample reads (these are pretrimmed, removing adapters via bbtools and primers via cutadapt)
path = "/scratch/mac9jc/crypto_micro/trimmed_seqs/run1"
#path <- "/project/tumi/carey/crypto_microbiome/trimmed_seqs/run1"
fns <- list.files(path)
fns = grep("_PT2.fq.gz",fns,value=TRUE)

# remove reads that couldn't be assigned to samples in the demultiplexing process
remove = grep("Undetermined_",fns,value=TRUE)
fns = fns[!fns %in% remove]

### Load forward and reverse reads
fastqs <- fns[grepl(".fq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_1_PT2", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_2_PT2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# remove samples if they have fewer than 10,000 reads prior to filtration
# this is VERY SLOW, so don't rerun everytime
remove_list = list()
for (i in 1:length(sample.names)) {
  f_seq = getSequences(fnFs[i]); r_seq = getSequences(fnRs[i])
  if (length(f_seq) < 10000 | length(r_seq) < 10000 ) {
   print(sample.names[i])
   remove_list = list(remove_list, i) }}
remove_list = unlist(remove_list)
# nto sure if this shortcut is working
#remove_list = which(c("130491008","130771018","131361012","neg1","neg2","neg3") %in% sample.names)
print(remove_list)
print('removing these samples from consideration because they have less than 10,000 reads prior to filtering')
print(sample.names[remove_list])
fnFs = fnFs[-remove_list]
print(length(fnFs))
fnRs = fnRs[-remove_list]
sample.names = sample.names[-remove_list]

print('about to make plot')
# visualize quality profile
plotQualityProfile(fnRs, aggregate = T)
ggsave("/project/tumi/carey/crypto_micro/plots/run1_reverse_reads_quality.pdf")
plotQualityProfile(fnFs, aggregate = T)
ggsave("/project/tumi/carey/crypto_micro/plots/run1_forward_reads_quality.pdf")

print('made plots')

### Filtering and trimming
# Make directory and filenames for the filtered fastqs
filt_path = "/project/tumi/carey/crypto_microbiome/filtered_seqs/run1"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fq.gz"))
# Filter
out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
              maxN=0, maxEE=c(3,6), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
saveRDS(out1, "/project/tumi/carey/crypto_microbiome/results/run1_out.rds")

print('made out file')

# In the above function, trimleft describes the number of bases to skip at the beginning of each read as (forward, reverse).
# truncLen is the same, except describes the base number to truncate at. Choice for truncLen is justified by the point at which
# error rates tend to start increasing.

filtpath <- filt_path # CHANGE ME to the directory containing your filtered fastq files
fns1 <- list.files(filtpath, full.names = TRUE)

# get forward read files for run 1
fnsF1 <- fns1[lapply('_F_',grepl,fns1)[[1]]]
# get reverse read files for run 1
fnsR1 <- fns1[lapply('_R_',grepl,fns1)[[1]]]

filtsF1 <- fnsF1[grepl("fq.gz$", fnsF1)] # CHANGE if different file extensions
filtsR1 <- fnsR1[grepl("fq.gz$", fnsR1)] # CHANGE if different file extensions

sample.names.F1 <- sapply(strsplit(basename(filtsF1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.names.R1 <- sapply(strsplit(basename(filtsR1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
names(filtsF1) <- sample.names.F1
names(filtsR1) <- sample.names.R1

print('start learning')
# learn error rates for each run, forward and reverse reads
set.seed(100)
errF1 <- learnErrors(filtsF1, nbases = 1e9, multithread=TRUE)
errR1 <- learnErrors(filtsR1, nbases = 1e9, multithread=TRUE)
save(errF1, errR1, file = "/project/tumi/carey/crypto_micro/err_for_plot_run1.RData")

p = plotErrors(errF1, nominalQ=TRUE)
ggsave("/project/tumi/carey/crypto_micro/plots/FER1.pdf",p)
p = plotErrors(errR1, nominalQ=TRUE)
ggsave("/project/tumi/carey/crypto_micro/plots/RER1.pdf",p)

print('start inference')

# sample inference on run 1
dadaFs1 <- dada(filtsF1, err=errF1, multithread=TRUE)
saveRDS(dadaFs1, "/project/tumi/carey/crypto_microbiome/results/run1_dadaFs1.rds")
dadaRs1 <- dada(filtsR1, err=errR1, multithread=TRUE)
saveRDS(dadaRs1, "/project/tumi/carey/crypto_microbiome/results/run1_dadaRs1.rds")
mergers1 <- mergePairs(dadaFs1, filtsF1, dadaRs1, filtsR1, maxMismatch = 1, verbose=TRUE) # merge forward and reverse reads
saveRDS(mergers1, "/project/tumi/carey/crypto_microbiome/results/run1_mergers.rds")
seqtab1 <- makeSequenceTable(mergers1) # construct sequence table
table(nchar(getSequences(seqtab1))) # Inspect distribution of sequence lengths
# remove chimeras
seqtab1 <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab1, "/project/tumi/carey/crypto_microbiome/results/run1_withtrimming_seqtab.rds")

print('these were in run 1 out file but not in mergers file')
temp = matrix(unlist((strsplit(rownames(out1), '_', fixed=FALSE))), ncol=6, byrow=TRUE)
print(setdiff(temp[,1], names(mergers1)))
# remove them
#l = list()
#for (x in setdiff(temp[,1], names(mergers1))) { l = list(l,which(temp[,1] == x))}
#l = unlist(l)
#out1 = out1[-l,]

print('start track reads')
# track reads through pipeline
getN <- function(x) sum(getUniques(x))
track1 <- cbind(out1, sapply(dadaFs1, getN), sapply(dadaRs1, getN), sapply(mergers1, getN), rowSums(seqtab1))
rownames(track1) = sample.names
colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
write.csv(track1, "/project/tumi/carey/crypto_microbiome/results/track_reads_run1.csv")
