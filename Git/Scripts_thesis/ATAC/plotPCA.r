library(DiffBind)

#BiocManager::install("profileplyr")
library(profileplyr)
# rm(list = ls())
# 
#load("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/data/minkangdiffbind_atac_h301_aml.RData")

pdf("pca.pdf", height = 7, width = 7)
p <- dba.plotPCA(diff, attributes = DBA_TREATMENT)
print(p)
dev.off()

pdf("pca2.pdf", height = 7, width = 7)
p <- dba.plotPCA(diff, attributes = DBA_TISSUE)
print(p)
dev.off()

pdf("correlation_heatmap.pdf", height = 7, width = 7)
p <- dba.plotHeatmap(diff)
print(p)
dev.off()

#only when degs are found 
pdf("venn_contrast_1.pdf", height = 7, width = 7 )
p <- dba.plotVenn(diff, contrast = 5, method = DBA_ALL_METHODS, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
print(p)
dev.off()

pdf("MA_c5.pdf", height = 7, width = 7)
p <- dba.plotMA(diff, contrast = 5)
print(p)
dev.off()

pdf("Vol_c5.pdf", height = 7, width = 7)
p <- dba.plotVolcano(diff, contrast = 5)
print(p)
dev.off()


pdf("box_pval.pdf", height = 7, width = 7)
sum(diff.DB$Fold<0)
sum(diff.DB$Fold>0)
p <- dba.plotBox(diff)
print(p)
dev.off()

# pdf("profiles_5.pdf", height = 7, width = 7)
# p <- DiffBind::dba.plotProfile(diff, merge = DBA_CONDITION)
# print(p)
# dev.off()