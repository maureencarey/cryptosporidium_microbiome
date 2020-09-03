#install.packages("BiocManager", repos='http://cran.us.r-project.org') # WORKS
#BiocManager::install("ShortRead", dependencies = TRUE)

library("dada2")
library("ShortRead")
library("phyloseq")
library("ggplot2")
sink('/project/tumi/carey/crypto_micro/dada2_R_outfile_r2.txt')
packageVersion('dada2')

# get raw sample reads (these are pretrimmed, removing adapters via bbtools and primers via cutadapt)
path2 = "/scratch/mac9jc/crypto_micro/trimmed_seqs/run2"
#path2 <- "/project/tumi/carey/crypto_microbiome/trimmed_seqs/run2"
fns2 <- list.files(path2)
fns2 = grep("_PT2.fq.gz",fns2,value=TRUE)

# remove reads that couldn't be assigned to samples in the demultiplexing process
remove = grep("Undetermined_",fns2,value=TRUE)
fns2 = fns2[!fns2 %in% remove]

### Load forward and reverse reads for datapool 2
fastqs2 <- fns2[grepl(".fq.gz$", fns2)]
fastqs2 <- sort(fastqs2) # Sort ensures forward/reverse reads are in same order
fnFs2 <- fastqs2[grepl("_1_PT2", fastqs2)] # Just the forward read files
fnRs2 <- fastqs2[grepl("__2_PT2", fastqs2)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names2 <- sapply(strsplit(fnFs2, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs2 <- file.path(path2, fnFs2)
fnRs2 <- file.path(path2, fnRs2)

# remove samples if they have fewer than 10,000 reads prior to filtration
# this is VERY SLOW, so don't rerun everytime
remove_list = list()
for (i in 1:length(sample.names2)) {
  f_seq = getSequences(fnFs2[i]); r_seq = getSequences(fnRs2[i])
  if (length(f_seq) < 10000 | length(r_seq) < 10000 ) {
   remove_list = list(remove_list, i) }}
remove_list = unlist(remove_list)
# replace with run 2 results and check if this line works in the other run
# remove_list = which(c("130491008","130771018","131361012","neg1","neg2","neg3") %in% sample.names2))
print('removing these samples from consideration because they have less than 10,000 reads prior to filtering')
print(sample.names2[remove_list])
fnFs2 = fnFs2[-remove_list]
print(length(fnFs2))
fnRs2 = fnRs2[-remove_list]
sample.names2 = sample.names2[-remove_list]

# visualize quality profile
plotQualityProfile(fnRs2, aggregate = T)
ggsave("/project/tumi/carey/crypto_micro/plots/run2_reverse_reads_quality.pdf")
plotQualityProfile(fnFs2, aggregate = T)
ggsave("/project/tumi/carey/crypto_micro/plots/run2_forward_reads_quality.pdf")

### Filtering and trimming
# Make directory and filenames for the filtered fastqs
filt_path2 = "/project/tumi/carey/crypto_microbiome/filtered_seqs/run2"
if(!file_test("-d", filt_path2)) dir.create(filt_path2)
filtFs2 <- file.path(filt_path2, paste0(sample.names2, "_F_filt.fq.gz"))
filtRs2 <- file.path(filt_path2, paste0(sample.names2, "_R_filt.fq.gz"))
# Filter
out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=c(250,200),
              maxN=0, maxEE=c(3,6), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
saveRDS(out2, "/project/tumi/carey/crypto_microbiome/results/run2_out.rds")

# In the above function, trimleft describes the number of bases to skip at the beginning of each read as (forward, reverse).
# truncLen is the same, except describes the base number to truncate at. Choice for truncLen is justified by the point at which
# error rates tend to start increasing.

filtpath2 <- filt_path2 # CHANGE ME to the directory containing your filtered fastq files
fns2 <- list.files(filtpath2, full.names = TRUE)

# get forward read files for run 1
fnsF2 <- fns2[lapply('_F_',grepl,fns2)[[1]]]
# get reverse read files for run 1
fnsR2 <- fns2[lapply('_R_',grepl,fns2)[[1]]]

filtsF2 <- fnsF2[grepl("fq.gz$", fnsF2)] # CHANGE if different file extensions
filtsR2 <- fnsR2[grepl("fq.gz$", fnsR2)] # CHANGE if different file extensions

sample.names.F2 <- sapply(strsplit(basename(filtsF2), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.names.R2 <- sapply(strsplit(basename(filtsR2), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
names(filtsF2) <- sample.names.F2
names(filtsR2) <- sample.names.R2

print('start learning')
# learn error rates for each run, forward and reverse reads
set.seed(100)
errF2 <- learnErrors(filtsF2, nbases = 1e9, multithread=TRUE)
errR2 <- learnErrors(filtsR2, nbases = 1e9, multithread=TRUE)
save(errF2, errR2, file = "/project/tumi/carey/crypto_micro/err_for_plot_r2.RData")

p = plotErrors(errF2, nominalQ=TRUE)
ggsave("/project/tumi/carey/crypto_micro/plots/FER2.pdf",p)
p = plotErrors(errR2, nominalQ=TRUE)
ggsave("/project/tumi/carey/crypto_micro/plots/RER2.pdf",p)

print('start inference')

# sample inference on run 2
dadaFs2 <- dada(filtsF2, err=errF2, multithread=TRUE)
saveRDS(dadaFs2, "/project/tumi/carey/crypto_microbiome/results/run2_dadaFs2.rds")
dadaRs2 <- dada(filtsR2, err=errR2, multithread=TRUE)
saveRDS(dadaRs2, "/project/tumi/carey/crypto_microbiome/results/run2_dadaRs2.rds")
mergers2 <- mergePairs(dadaFs2, filtsF2, dadaRs2, filtsR2, maxMismatch = 1, verbose=TRUE) # merge forward and reverse reads
saveRDS(mergers2, "/project/tumi/carey/crypto_microbiome/results/run2_mergers.rds")
seqtab2 <- makeSequenceTable(mergers2) # construct sequence table
table(nchar(getSequences(seqtab2))) # Inspect distribution of sequence lengths
# remove chimeras
seqtab2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab2, "/project/tumi/carey/crypto_microbiome/results/run2_withtrimming_seqtab.rds")

print('these were in run 2 out file but not in mergers file')
temp = matrix(unlist((strsplit(rownames(out2), '_', fixed=FALSE))), ncol=6, byrow=TRUE)
print(setdiff(temp[,1], names(mergers2)))
# remove them
#l = list()
#for (x in setdiff(temp[,1], names(mergers2))) { l = list(l,which(temp[,1] == x))}
#l = unlist(l)
#out2 = out2[-l,]

print('start track reads')
# track reads through pipeline
getN <- function(x) sum(getUniques(x))
track2 <- cbind(out2, sapply(dadaFs2, getN), sapply(dadaRs2, getN), sapply(mergers2, getN), rowSums(seqtab2))
rownames(track2) = sample.names2
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
write.csv(track2, "/project/tumi/carey/crypto_microbiome/results/track_reads_run2.csv")

