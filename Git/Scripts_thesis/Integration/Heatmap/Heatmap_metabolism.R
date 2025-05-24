# heatmap metbolism 

library(pheatmap)
library(tidyverse)

#load data

t_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_trans.csv")
rownames(t_z_score) <- t_z_score$X  
t_z_score$X <- NULL

p_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_proteomics.csv")
rownames(p_z_score) <- p_z_score$X  
p_z_score$X <- NULL

# gene list metabolism
glycolysis_genes <- c("ACADSB", "ACSF3", "ALDH3A2", "CPT1A", 
                      "PCCA", "SLC27A1", "SLC1A5", "SLC7A1", "CYCS", 
                       "FDFT1", "ATP5ME","ATP5MF", "ATP5F1A", "ALDOA", "ENO1", "GAPDH", "GPI", "HK1", "HK2", "LDHA", 
                      "PFKM", "PGAM1", "PKM", "PKMYT1", "SLC2A1", "ACSL4", "DHCR24", "FADS2", "FASN", 
                       "FDPS", "HMGCS1", "IDI1", 
                      "MVD", "MVK", "SCD", "SOAT1", "SQLE", "SREBF1", "ACACA")


glycolysis_genes <- sort(glycolysis_genes)

# filter data on gene list 
transcriptomics_filtered <- t_z_score %>%
  filter(row.names(t_z_score) %in% glycolysis_genes)

proteomics_filtered <- p_z_score %>%
  filter(row.names(p_z_score) %in% glycolysis_genes)

print(colnames(proteomics_filtered))
# change condition names and orders 
selected_conditions_t <- c("aCD3", "aCD3.1", "aCD3.2", "aCD3.3", "aCD3.4",'aCD3_aCD27', "aCD3_aCD27.1", "aCD3_aCD27.2", "aCD3_aCD27.3", "aCD3_aCD27.4", "aCD3_a4_1BB", "aCD3_a4_1BB.1", "aCD3_a4_1BB.2", "aCD3_a4_1BB.3", "aCD3_a4_1BB.4", "aCD3_aCD28", "aCD3_aCD28.1", "aCD3_aCD28.2", "aCD3_aCD28.3", "aCD3_aCD28.4", "aCD3_aCD27_aCD28", "aCD3_aCD27_aCD28.1", "aCD3_aCD27_aCD28.2", "aCD3_aCD27_aCD28.3", "aCD3_aCD27_aCD28.4", "aCD3_aCD28_a4_1BB", "aCD3_aCD28_a4_1BB.1", "aCD3_aCD28_a4_1BB.2", "aCD3_aCD28_a4_1BB.3", "aCD3_aCD28_a4_1BB.4") # Vul hier je gewenste condities in
selected_conditions_p <- c("X2.1.aCD3", "X2.2.aCD3", "X2.3.aCD3", "X3.1.aCD3.aCD27", "X3.2.aCD3.aCD27", "X3.3.aCD3.aCD27", "X4.1.aCD3.a4_1BB", "X4.2.aCD3.a4_1BB", "X4.3.aCD3.a4_1BB", "X5.1.aCD3.aCD28", "X5.2.aCD3.aCD28", "X5.3.aCD3.aCD28", "X6.1.aCD3.aCD27.aCD28", "X6.2.aCD3.aCD27.aCD28", "X6.3.aCD3.aCD27.aCD28", "X7.1.aCD3.aCD28.a4_1BB", "X7.2.aCD3.aCD28.a4_1BB", "X7.3.aCD3.aCD28.a4_1BB")

# order 
transcriptomics_filtered <- transcriptomics_filtered[glycolysis_genes, ,drop = FALSE]
proteomics_filtered <- proteomics_filtered[glycolysis_genes, ,drop = FALSE]

#filter row 
transcriptomics_filtered <- transcriptomics_filtered[, selected_conditions_t, drop = FALSE]
proteomics_filtered <- proteomics_filtered[, selected_conditions_p, drop = FALSE]

# choose colours 
heatmap_colors_t <- colorRampPalette(c("blue", "black", "yellow"))(50)
heatmap_colors_p <- colorRampPalette(c("blue", "white", "red"))(50)




# heatmap transcriptomic
p1 <- pheatmap(transcriptomics_filtered, 
               cluster_rows = FALSE, cluster_cols = FALSE,  # Clustert genen maar behoudt de condities volgorde
               scale = "row",  
               main = "Transcriptomics Heatmap - Metabolism", 
               color = heatmap_colors_t)

# heatmap proteomic 
p2 <- pheatmap(proteomics_filtered, 
               cluster_rows = FALSE, cluster_cols = FALSE,  
               scale = "row", 
               main = "Proteomics Heatmap - Metabolism", 
               color = heatmap_colors_p)

# put heatmap next to eachother 
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

heatmap1 <- pheatmap(transcriptomics_filtered, silent = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, main = "Transcriptomic heatmap (z-scores) - Metabolism", color = heatmap_colors_t, fontsize = 7, fontsize_row = 10)$gtable
heatmap2 <- pheatmap(proteomics_filtered, silent = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, main = "Proteomic heatmap (z-scores) - Metabolism", color = heatmap_colors_p, fontsize = 7, fontsize_row = 10)$gtable

p3 <- gridExtra::grid.arrange(heatmap1, heatmap2, ncol = 2)

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/heatmaps_genes_felix/new_metabol/"
output_base <- paste0(output_dir, "Heatmap_metabol_trans")

# Save
ggsave(paste0(output_base, ".png"), plot = p1, width = 8, height = 6, dpi = 300)
ggsave(paste0(output_base, ".svg"), plot = p1, width = 8, height = 6)
ggsave(paste0(output_base, ".eps"), plot = p1, width = 8, height = 6, device = cairo_ps)



output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/heatmaps_genes_felix/new_metabol/"
output_base <- paste0(output_dir, "Heatmap_metabol_prot")

# save
ggsave(paste0(output_base, ".png"), plot = p2, width = 8, height = 6, dpi = 300)

ggsave(paste0(output_base, ".svg"), plot = p2, width = 8, height = 6)
ggsave(paste0(output_base, ".eps"), plot = p2, width = 8, height = 6, device = cairo_ps)


output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/heatmaps_genes_felix/new_metabol/"
output_base <- paste0(output_dir, "Heatmap_metabol_both")

# Save
ggsave(paste0(output_base, ".png"), plot = p3, width = 8, height = 6, dpi = 300)
ggsave(paste0(output_base, ".svg"), plot = p3, width = 8, height = 6)

ggsave(paste0(output_base, ".eps"), plot = p3, width = 8, height = 6, device = cairo_ps)

