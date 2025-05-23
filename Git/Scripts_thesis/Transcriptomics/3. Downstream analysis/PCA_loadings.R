library(DESeq2)

# Perform variance-stabilizing transformation (VST)
vsd <- vst(dds_csp, blind = FALSE)

# Extract the transformed data
vsd_datacsp <- assay(vsd)  

# Perform PCA manually on the transposed data (samples as rows, genes as columns)
pca_rescsp <- prcomp(t(vsd_datacsp), scale. = TRUE)

# Extract the loadings (gene contributions to the PCs)
loadings <- pca_rescsp$rotation  

# Identify top contributing genes for PC1 and PC2 (Top 10 by absolute loading values)
top_genes_pc1 <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:10]
top_genes_pc2 <- sort(abs(loadings[, 2]), decreasing = TRUE)[1:10]

# Calculate percentage contributions for these top genes
top_genes_pc1_contributions <- (top_genes_pc1^2 / sum(top_genes_pc1^2)) * 100
top_genes_pc2_contributions <- (top_genes_pc2^2 / sum(top_genes_pc2^2)) * 100

# Prepare a data frame with gene names and their percentage contribution to PC1 and PC2
top_genes_pc1_names <- names(top_genes_pc1)
top_genes_pc2_names <- names(top_genes_pc2)

top_genes_df <- data.frame(
  Gene = c(top_genes_pc1_names, top_genes_pc2_names),
  Contribution = c(top_genes_pc1_contributions, top_genes_pc2_contributions),
  PC = rep(c("PC1", "PC2"), each = 10)
)

print(top_genes_df)

#contribution all genes manually check
loadings_pc1 <- pca_rescsp$rotation[,1]
contribution_pc1 <- loadings_pc1^2/sum(loadings_pc1^2)*100
contribution_pc1["ENSG00000108424"]