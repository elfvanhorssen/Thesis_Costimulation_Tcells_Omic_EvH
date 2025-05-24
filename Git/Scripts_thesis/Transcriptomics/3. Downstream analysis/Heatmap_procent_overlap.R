library(pheatmap)

# get the names of the comparisons 
comparisons <- names(deseq_res)
n <- length(comparisons)

# Make an empthy matrix for the overlap % 
overlap_matrix <- matrix(0, nrow = n, ncol = n)
colnames(overlap_matrix) <- comparisons
rownames(overlap_matrix) <- comparisons

for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      res_i <- deseq_res[[comparisons[i]]]
      res_j <- deseq_res[[comparisons[j]]]
      
      # get gene ids and make tresholds 
      genes_up_i <- rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 & res_i$log2FoldChange > 1 ]
      genes_up_j <- rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 & res_j$log2FoldChange > 1 ]
      
      genes_down_i <- rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 & res_i$log2FoldChange < -1 ]
      genes_down_j <- rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 & res_j$log2FoldChange < -1 ]
      
      # count common up and downregulated genes 
      common_up <- sum(genes_up_i %in% genes_up_j)
      common_down <- sum(genes_down_i %in% genes_down_j)
      
      # get the total deg or deps per comparison 
      total_genes <- length(unique(c(
        rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 ],
        rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 ]
      )))
      
      # calculate overlap. 
      if(total_genes > 0) {
        overlap_matrix[i, j] <- (common_up + common_down) / total_genes * 100
      } else {
        overlap_matrix[i, j] <- 0
      }
    }
  }
}

# diagonal is 100% 
diag(overlap_matrix) <- 100

# Define comparisons 
selected_comparisons <- c("aCD3_aCD27_vs_aCD3", 
                          "aCD3_aCD28_vs_aCD3", 
                          "aCD3_a4_1BB_vs_aCD3", 
                          "aCD3_aCD27_aCD28_vs_aCD3", 
                          "aCD3_aCD28_a4_1BB_vs_aCD3")


# Filter matrix on the chosen comparisons 
overlap_matrix_filtered <- overlap_matrix[selected_comparisons, selected_comparisons]

#order matrix 
desired_order <- c("aCD3_aCD27_vs_aCD3", 
                   "aCD3_a4_1BB_vs_aCD3", 
                   "aCD3_aCD28_vs_aCD3",
                   "aCD3_aCD27_aCD28_vs_aCD3", 
                   "aCD3_aCD28_a4_1BB_vs_aCD3")

# rows and columns in the same order 
overlap_matrix_filtered <- overlap_matrix_filtered[desired_order, desired_order]

#make vector with breaks 
my_breaks <- seq(0, 100, by = 1)

# choose colors 
my_colors <- colorRampPalette(c("white", "#FFC273"))(length(my_breaks) - 1)

#heatmap overlap
p <- pheatmap(
  overlap_matrix_filtered,
  color            = my_colors,
  breaks           = my_breaks,
  cluster_rows     = FALSE,     
  cluster_cols     = FALSE,
  display_numbers  = TRUE,      
  number_format    = "%.1f",    
  fontsize_number  = 12,        
  number_color     = "black",   
  fontsize         = 12,        
  fontsize_row     = 12,        
  fontsize_col     = 12,        
  angle_col        = 45,        
  border_color     = "grey90",  
  main             = "Overlap (%) of DE genes between Comparisons compared to aCD3"
)

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/procent_heatmap_trans"
output_base <- paste0(output_dir, "procent_DE_transcriptomics")

# Save
ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)

ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)

ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)


