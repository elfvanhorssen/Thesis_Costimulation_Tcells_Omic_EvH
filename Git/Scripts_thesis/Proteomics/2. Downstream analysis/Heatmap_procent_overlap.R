result_list_log <- readRDS("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/result_list_log_prot.rds")

library(pheatmap)
# get the names for condition out of limma 
comparisons <- names(result_list_log)
n <- length(comparisons)

# Make an empthy matrix for the overlap % 
overlap_matrix <- matrix(0, nrow = n, ncol = n)
colnames(overlap_matrix) <- comparisons
rownames(overlap_matrix) <- comparisons

# choose threshold 
pval_threshold <- 0.05
logfc_threshold <- 0.5

# Loop trough every comparison 
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      res_i <- result_list_log[[comparisons[i]]]
      res_j <- result_list_log[[comparisons[j]]]
      
      # get gene/proteinID
      genes_up_i <- rownames(res_i)[ !is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC > logfc_threshold ]
      genes_up_j <- rownames(res_j)[ !is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC > logfc_threshold ]
      
      genes_down_i <- rownames(res_i)[ !is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC < -logfc_threshold ]
      genes_down_j <- rownames(res_j)[ !is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC < -logfc_threshold ]
      
      # count common up and downregulated genes 
      common_up <- sum(genes_up_i %in% genes_up_j)
      common_down <- sum(genes_down_i %in% genes_down_j)
      
      # get the total deg or deps per comparison 
      total_genes <- length(unique(c(
        rownames(res_i)[ !is.na(res_i$P.Value) & res_i$P.Value < pval_threshold ],
        rownames(res_j)[ !is.na(res_j$P.Value) & res_j$P.Value < pval_threshold ]
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
selected_comparisons <- c("aCD3_aCD27 - aCD3", 
                          
                          "aCD3_aCD28 - aCD3",
                          "aCD3_a4_1BB - aCD3", 
                          "aCD3_aCD27_aCD28 - aCD3", 
                          "aCD3_aCD28_a4_1BB - aCD3")

# Filter matrix on the chosen comparisons 
overlap_matrix_filtered <- overlap_matrix[selected_comparisons, selected_comparisons]

# order
desired_order <- c("aCD3_aCD27 - aCD3", 
                   "aCD3_a4_1BB - aCD3", 
                   "aCD3_aCD28 - aCD3",
                   
                   "aCD3_aCD27_aCD28 - aCD3", 
                   "aCD3_aCD28_a4_1BB - aCD3")
overlap_matrix_filtered <- overlap_matrix_filtered[desired_order, desired_order]

#make an vector with breaks 
my_breaks <- seq(0, 100, by = 1)

# colors 
my_colors <- colorRampPalette(c("white", "orange"))(length(my_breaks) - 1)

#make heatmap without clustering 
p <- pheatmap(
  overlap_matrix_filtered,
  color            = my_colors,
  breaks           = my_breaks,
  cluster_rows     = FALSE,     
  cluster_cols     = FALSE,
  display_numbers  = TRUE,      # show overlap in % 
  number_format    = "%.1f",    # 1 decimal 
  fontsize_number  = 12,        
  number_color     = "black",   
  fontsize         = 12,        
  fontsize_row     = 12,       
  fontsize_col     = 12,        
  angle_col        = 45,        
  border_color     = "grey90",  
  main             = "Overlap (%) of DE proteins between Comparisons compared to aCD3"
)

p

#save 
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/procent_heatmap_prot/"
output_base <- paste0(output_dir, "procent_DE_proteomics")


ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)

ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)

ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)
