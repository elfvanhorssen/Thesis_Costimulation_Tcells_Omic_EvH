result_list_log <- readRDS("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/result_list_log_prot.rds")
library(pheatmap)

# load comparison names and initialize matrix
comparisons <- names(result_list_log)
n <- length(comparisons)
overlap_matrix <- matrix(0, nrow = n, ncol = n,
                         dimnames = list(comparisons, comparisons))

# set thresholds
pval_threshold  <- 0.05
logfc_threshold <- 0.5

# calculate overlap (%) between each pair
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i != j) {
      res_i <- result_list_log[[comparisons[i]]]
      res_j <- result_list_log[[comparisons[j]]]
      
      genes_up_i   <- rownames(res_i)[!is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC >  logfc_threshold]
      genes_up_j   <- rownames(res_j)[!is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC >  logfc_threshold]
      genes_down_i <- rownames(res_i)[!is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC < -logfc_threshold]
      genes_down_j <- rownames(res_j)[!is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC < -logfc_threshold]
      
      common_up   <- sum(genes_up_i   %in% genes_up_j)
      common_down <- sum(genes_down_i %in% genes_down_j)
      
      total_genes <- length(unique(c(
        rownames(res_i)[!is.na(res_i$P.Value) & res_i$P.Value < pval_threshold],
        rownames(res_j)[!is.na(res_j$P.Value) & res_j$P.Value < pval_threshold]
      )))
      
      overlap_matrix[i, j] <- if (total_genes > 0) {
        (common_up + common_down) / total_genes * 100
      } else {
        0
      }
    }
  }
}

# set self-overlap to 100%
diag(overlap_matrix) <- 100

# filter and reorder for heatmap
selected <- c("aCD3_aCD27 - aCD3",
              "aCD3_aCD28 - aCD3",
              "aCD3_a4_1BB - aCD3",
              "aCD3_aCD27_aCD28 - aCD3",
              "aCD3_aCD28_a4_1BB - aCD3")
overlap_subset <- overlap_matrix[selected, selected]

# set breaks and color palette
breaks <- seq(0, 100, by = 1)
colors <- colorRampPalette(c("white", "orange"))(length(breaks) - 1)

# plot heatmap without clustering
p <- pheatmap(
  overlap_subset,
  color           = colors,
  breaks          = breaks,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = TRUE,
  number_format   = "%.1f",
  fontsize_number = 12,
  number_color    = "black",
  fontsize        = 12,
  angle_col       = 45,
  border_color    = "grey90",
  main            = "overlap (%) of DE proteins vs aCD3"
)

# save heatmap in multiple formats
out_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/procent_heatmap_prot/"
base    <- file.path(out_dir, "procent_DE_proteomics")
ggsave(paste0(base, ".png"), p, width = 8, height = 6, dpi = 300)
ggsave(paste0(base, ".svg"), p, width = 8, height = 6)
ggsave(paste0(base, ".eps"), p, width = 8, height = 6, device = cairo_ps)
