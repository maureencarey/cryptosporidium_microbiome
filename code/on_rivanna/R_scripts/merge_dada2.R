#install.packages("BiocManager", repos='http://cran.us.r-project.org') # WORKS
#BiocManager::install("ShortRead", dependencies = TRUE)

library("dada2")
library("ShortRead")
library("phyloseq")
sink('/project/tumi/carey/crypto_micro/dada2_R_outfile_merge.txt')
packageVersion('dada2')

# load both sequence tables, merge
st1 <- readRDS("/project/tumi/carey/crypto_microbiome/results/run1_withtrimming_seqtab.rds")
st2 <- readRDS("/project/tumi/carey/crypto_microbiome/results/run2_withtrimming_seqtab.rds")

# append _run1 and _run2 to each sample name
rownames(st1) <- paste0(rownames(st1),'_run1')
rownames(st2) <- paste0(rownames(st2),'_run2')

st <- mergeSequenceTables(st1, st2)
# remove sequences far outside the desired length
# st <- st[,nchar(colnames(st)) %in% seq(272,282)]

# collapse nomismatch sequences, which are identical other than being of differing lengths
st = collapseNoMismatch(st)

print('start assigning taxonomy')
# assign taxonomy
taxa <- assignTaxonomy(st,"/project/tumi/carey/crypto_microbiome/rdp_train_set_16.fa")

### transfer to phyloseq
#prep
samples.out <- rownames(st)
# load and format metadata
meta <- read.table("/project/tumi/carey/crypto_microbiome/meta_16_12_5.csv",header=TRUE, sep='\t')
# append _run1 and _run2 to sample id to match sequence table rows
a = paste0(meta$sample.id,'_run')
meta$sample.id <- paste0(a,meta$Library)
rownames(meta) <- meta$sample.id

#match rownames for merging
otu_table_use = otu_table(st, taxa_are_rows = FALSE)
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos1_run1'] = 'pos1_run1'
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos2_run1'] = 'pos2_run1'
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos3_run1'] = 'pos3_run1'
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos4_run1'] = 'pos4_run1'
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos5_run1'] = 'pos5_run1'
rownames(otu_table_use)[rownames(otu_table_use) == 'HM-782D-Pos6_run1'] = 'pos6_run1'

sample_data_use = sample_data(meta)
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos1)_run1'] = 'pos1_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos2)_run1'] = 'pos2_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos3)_run1'] = 'pos3_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos4)_run1'] = 'pos4_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos5)_run1'] = 'pos5_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos6)_run1'] = 'pos6_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg1)_run1'] = 'neg1_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg2)_run1'] = 'neg2_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg3)_run1'] = 'neg3_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg4)_run1'] = 'neg4_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg5)_run1'] = 'neg5_run1'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg6)_run1'] = 'neg6_run1'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos7)_run2'] = 'pos1_run2'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos8)_run2'] = 'pos2_run2'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos9)_run2'] = 'pos3_run2'
rownames(sample_data_use)[rownames(sample_data_use) ==' HM-782D (pos10)_run2'] = 'pos4_run2'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg7)_run2'] = 'neg1_run2'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg8)_run2'] = 'neg2_run2'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg9)_run2'] = 'neg3_run2'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg10)_run2'] = 'neg4_run2'
rownames(sample_data_use)[rownames(sample_data_use) =='NTC (neg11)_run2'] = 'neg5_run2'

# drop these samples because they were pruned due to having too few reads
# DOUBLE CHECK IF RERUNNING - these samples will be printed in in dada2_R_outfile_r1.txt ~line3
drop_these = c('130491008_run1','130771018_run1','131361012_run1','neg1_run1','neg2_run1','neg3_run1',
               'neg4_run1','neg5_run1','neg6_run1')
sample_data_use = sample_data_use[!(rownames(sample_data_use) %in% drop_these),]

# drop these samples because they were pruned due to having too few reads 
# DOUBLE CHECK IF RERUNNING - these samples will be printed in in dada2_R_outfile_r2.txt ~line3
drop_these = c('130751006_run2','131971007_run2','132811005_run2','neg1_run2','neg2_run2','neg3_run2','neg4_run2','neg5_run2')
sample_data_use = sample_data_use[!(rownames(sample_data_use) %in% drop_these),]

print('nothing should print in these lists:')
print(setdiff(rownames(sample_data_use),rownames(otu_table_use)))
print(setdiff(rownames(otu_table_use),rownames(sample_data_use)))

#phyloseq
ps <- phyloseq(otu_table_use,sample_data_use,tax_table(taxa))
# save phyloseq object to load later
saveRDS(ps, "/project/tumi/carey/crypto_microbiome/results/phyloseq_withtrim_withcollapse_object.rds")

# # DESeq differential abundance comparisons
# # Global site comparison first
# library(DESeq2)
# sample_data(ps1)$Clinical.type = as.factor(sample_data(ps1)$Clinical.type)
# deseq2_obj <- phyloseq_to_deseq2(ps1,~Clinical.type)
# #sizefactors <- estimateSizeFactors(deseq2_obj,type='iterate')
# 
# #sizefactors <- estimateSizeFactors(deseq2_input,type='iterate')
# deseq2_output <- DESeq(deseq2_obj, test="Wald", fitType="parametric")
# deseq2_output <- DESeqDataSet(deseq2_input
# deseq2_results <- results(deseq2_output, cooksCutoff = FALSE)
# 
# alpha = 0.05
# sigtab = deseq2_results[which(deseq2_results$padj < alpha), ]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps1)[rownames(sigtab), ], "matrix"))
# head(sigtab_d13)
