result_list_log <- readRDS("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/result_list_log_prot.rds")

library(pheatmap)

# -------------------------------
# 1. Overlap tussen comparisons berekenen (limma)
# -------------------------------

# Verkrijg de namen van de comparisons uit je limma-resultaten
comparisons <- names(result_list_log)
n <- length(comparisons)

# Maak een lege matrix voor de overlap (%) tussen comparisons
overlap_matrix <- matrix(0, nrow = n, ncol = n)
colnames(overlap_matrix) <- comparisons
rownames(overlap_matrix) <- comparisons

# Stel drempelwaarden in (pas deze eventueel aan)
pval_threshold <- 0.05
logfc_threshold <- 0.5

# Loop door alle combinations van comparisons
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      res_i <- result_list_log[[comparisons[i]]]
      res_j <- result_list_log[[comparisons[j]]]
      
      # Haal de gene/proteïne IDs op voor de criteria: 
      # Let op: bij limma gebruiken we vaak 'adj.P.Val' en 'logFC'
      genes_up_i <- rownames(res_i)[ !is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC > logfc_threshold ]
      genes_up_j <- rownames(res_j)[ !is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC > logfc_threshold ]
      
      genes_down_i <- rownames(res_i)[ !is.na(res_i$adj.P.Val) & res_i$P.Value < pval_threshold & res_i$logFC < -logfc_threshold ]
      genes_down_j <- rownames(res_j)[ !is.na(res_j$adj.P.Val) & res_j$P.Value < pval_threshold & res_j$logFC < -logfc_threshold ]
      
      # Tel het aantal gemeenschappelijke Up- en Downregulated genes
      common_up <- sum(genes_up_i %in% genes_up_j)
      common_down <- sum(genes_down_i %in% genes_down_j)
      
      # Bepaal het totaal aantal unieke significante genes in beide comparisons
      total_genes <- length(unique(c(
        rownames(res_i)[ !is.na(res_i$P.Value) & res_i$P.Value < pval_threshold ],
        rownames(res_j)[ !is.na(res_j$P.Value) & res_j$P.Value < pval_threshold ]
      )))
      
      # Bereken de overlap in percentage (voorkom deling door 0)
      if(total_genes > 0) {
        overlap_matrix[i, j] <- (common_up + common_down) / total_genes * 100
      } else {
        overlap_matrix[i, j] <- 0
      }
    }
  }
}

# Zet de diagonaal op 100% (elke comparison overlapt volledig met zichzelf)
diag(overlap_matrix) <- 100

# -------------------------------
# 2. Filteren en herordenen voor de heatmap
# -------------------------------

# Definieer de gewenste comparisons (pas deze aan naar jouw situatie)
selected_comparisons <- c("aCD3_aCD27 - aCD3", 
                          
                          "aCD3_aCD28 - aCD3",
                          "aCD3_a4_1BB - aCD3", 
                          "aCD3_aCD27_aCD28 - aCD3", 
                          "aCD3_aCD28_a4_1BB - aCD3")

# Filter de overlap matrix zodat alleen de geselecteerde comparisons overblijven
overlap_matrix_filtered <- overlap_matrix[selected_comparisons, selected_comparisons]

# Herordenen van de matrix volgens de gewenste volgorde
desired_order <- c("aCD3_aCD27 - aCD3", 
                   "aCD3_a4_1BB - aCD3", 
                   "aCD3_aCD28 - aCD3",
                   
                   "aCD3_aCD27_aCD28 - aCD3", 
                   "aCD3_aCD28_a4_1BB - aCD3")
overlap_matrix_filtered <- overlap_matrix_filtered[desired_order, desired_order]

# -------------------------------
# 3. Kleur- en schaalinstellingen
# -------------------------------
# Maak een vector van breaks van 0 t/m 100
my_breaks <- seq(0, 100, by = 1)

# Genereer een kleurenpalet (van wit naar coral, pas eventueel aan)
my_colors <- colorRampPalette(c("white", "orange"))(length(my_breaks) - 1)

# -------------------------------
# 4. Heatmap plotten (zonder clustering)
# -------------------------------
p <- pheatmap(
  overlap_matrix_filtered,
  color            = my_colors,
  breaks           = my_breaks,
  cluster_rows     = FALSE,     
  cluster_cols     = FALSE,
  display_numbers  = TRUE,      # Toon de overlap (%) in de cellen
  number_format    = "%.1f",    # Afronden op 1 decimaal
  fontsize_number  = 12,        # Tekstgrootte voor de getallen
  number_color     = "black",   # Kleur van de getallen
  fontsize         = 12,        # Algemene tekstgrootte
  fontsize_row     = 12,        # Rijlabels
  fontsize_col     = 12,        # Kolomlabels
  angle_col        = 45,        # Draai kolomlabels 45 graden
  border_color     = "grey90",  # Lichte rand om de cellen
  main             = "Overlap (%) of DE proteins between Comparisons compared to aCD3"
)

p
# -------------------------------
# 5. Opslaan van de heatmap
# -------------------------------
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/procent_heatmap_prot/"
output_base <- paste0(output_dir, "procent_DE_proteomics")

# Sla op als PNG (hoge resolutie)
ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)

# Sla op als SVG
ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)

# Sla op als EPS (voor vectorformaat)
ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)
