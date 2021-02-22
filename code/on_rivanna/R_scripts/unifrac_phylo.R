#install.packages("BiocManager", repos='http://cran.us.r-project.org') # WORKS
#BiocManager::install("dada2", dependencies = TRUE)

library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")
library("ggplot2")
sink('/project/tumi/carey/crypto_micro/phylo_out.txt')
packageVersion('dada2')

load("/home/mac9jc/R_obj.Rdata")

ps_rare = ps_rare_for_beta
sequences = dada2::getSequences(otu_table(ps_rare))
names(sequences) = sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
print('aligned')
dm <- dist.ml(phang.align)
print('dist ml')
treeNJ <- NJ(dm)
print('NJ')
fit = pml(treeNJ, data=phang.align)
print('pml')
fitGTR <- update(fit, k=4, inv=0.2)
print('update')
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
print('optim.pml')
print(fitGTR)
detach("package:DECIPHER", unload=TRUE)
detach("package:phangorn", unload=TRUE)
phy_tree(ps_rare) = phy_tree(fitGTR$tree) 
save(list = c("ps_rare"),file = paste0(paste0("/home/mac9jc/new_R_obj_" , format(Sys.time(), "%Y_%m_%d")),".RData"))

# calc distance and do ordination
dist_obj = phyloseq::distance(ps_rare, method = "wunifrac") #uunifrac for unweighted
print('dist ibj')
ord_obj = ordinate(ps_rare, method = "PCoA", distance = dist_obj) # default method is DCA
print('ord obj')

#plot_ordination(ps_rare, ord_obj, type = "samples")
#
p = plot_ordination(
  physeq = ps_rare,
  ordination = ord_obj,
  color = "Clinical.type",
  shape = "Number.of.crypto.positive.events"
)   + 
  scale_color_manual(values = c("#666666","#FF7F00","#1F78B4"),
                     labels = c("12" = "diarrheal","10" = "subclinical","1" = "pre-detection") )+ # 12 is diarrhea, 1,2 are pre
  theme(plot.margin=unit(c(.5,.5,.5,.5),"cm"))  + 
  labs(color = "infection", shape = "Number of infections")
ggsave("/home/mac9jc/Fig4A_pcoa_mirpur_unifrac.png", 
       plot=p, width = 5.4,height = 3,units = "in",dpi = 600)
