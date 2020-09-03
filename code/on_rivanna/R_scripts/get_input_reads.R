library("dada2")
library("ShortRead")
sink('/project/tumi/carey/crypto_micro/R_outfile_get_reads.txt')

# get raw sample reads (these are pretrimmed, removing adapters via bbtools and primers via cutadapt)
path1 = "/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool1/fastqs"
fns1 <- list.files(path1)
path2 = "/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool2/fastqs"
fns2 <- list.files(path2)
fns1 = grep(".fastq.gz",fns1,value=TRUE)
fns2 = grep(".fastq.gz",fns2,value=TRUE)

# remove reads that couldn't be assigned to samples in the demultiplexing process
remove = grep("Undetermined_",fns1,value=TRUE)
fns1 = fns1[!fns1 %in% remove]
remove = grep("Undetermined_",fns2,value=TRUE)
fns2 = fns2[!fns2 %in% remove]

### Load forward and reverse reads
fastqs1 <- fns1[grepl(".fastq.gz$", fns1)]
fastqs1 <- sort(fastqs1) # Sort ensures forward/reverse reads are in same order
fnFs1 <- fastqs1[grepl("_R1_001", fastqs1)] # Just the forward read files
fnRs1 <- fastqs1[grepl("_R2_001", fastqs1)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names1 <- sapply(strsplit(fnFs1, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs1 <- file.path(path1, fnFs1)
fnRs1 <- file.path(path1, fnRs1)

fastqs2 <- fns2[grepl(".fastq.gz$", fns2)]
fastqs2 <- sort(fastqs2) # Sort ensures forward/reverse reads are in same order
fnFs2 <- fastqs2[grepl("_R1_001", fastqs2)] # Just the forward read files
fnRs2 <- fastqs2[grepl("_R2_001", fastqs2)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names2 <- sapply(strsplit(fnFs2, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs2 <- file.path(path2, fnFs2)
fnRs2 <- file.path(path2, fnRs2)

df = data.frame(sample = sample.names1, forward = NA, reverse = NA)
df2 = data.frame(sample = sample.names2, forward = NA, reverse = NA)

for (i in 1:length(sample.names1)) {
  f_seq = getSequences(fnFs1[i]); r_seq = getSequences(fnRs1[i])
  df[i,2] = length(f_seq); df[i,3] = length(r_seq)}
for (i in 1:length(sample.names2)) {
  f_seq = getSequences(fnFs2[i]); r_seq = getSequences(fnRs2[i])
  df2[i,2] = length(f_seq); df2[i,3] = length(r_seq)}

write.csv(df, "/project/tumi/carey/crypto_micro/starting_seqs_run1.csv")
write.csv(df2, "/project/tumi/carey/crypto_micro/starting_seqs_run2.csv")
