library(DiffBind)
# rm(list = ls())
# 
# load("diffbind_atac_h301_aml.RData")
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


