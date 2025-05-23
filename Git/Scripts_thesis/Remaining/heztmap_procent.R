library(pheatmap)

# -------------------------------
# 1. Overlap tussen comparisons berekenen
# -------------------------------

# Verkrijg de namen van de comparisons
comparisons <- names(deseq_res)
n <- length(comparisons)

# Maak een lege matrix voor de overlap (%)
overlap_matrix <- matrix(0, nrow = n, ncol = n)
colnames(overlap_matrix) <- comparisons
rownames(overlap_matrix) <- comparisons

for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      res_i <- deseq_res[[comparisons[i]]]
      res_j <- deseq_res[[comparisons[j]]]
      
      # Haal de gene IDs op (gebruik rownames) voor de criteria:
      genes_up_i <- rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 & res_i$log2FoldChange > 1 ]
      genes_up_j <- rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 & res_j$log2FoldChange > 1 ]
      
      genes_down_i <- rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 & res_i$log2FoldChange < -1 ]
      genes_down_j <- rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 & res_j$log2FoldChange < -1 ]
      
      # Tel het aantal gemeenschappelijke Up- en Downregulated genen
      common_up <- sum(genes_up_i %in% genes_up_j)
      common_down <- sum(genes_down_i %in% genes_down_j)
      
      # Bepaal het totaal aantal unieke significante genen in beide comparisons
      total_genes <- length(unique(c(
        rownames(res_i)[ !is.na(res_i$padj) & res_i$padj < 0.01 ],
        rownames(res_j)[ !is.na(res_j$padj) & res_j$padj < 0.01 ]
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

# Zet de diagonaal op 100% (een comparison overlapt volledig met zichzelf)
diag(overlap_matrix) <- 100

# -------------------------------
# 2. Filteren en herordenen voor de heatmap
# -------------------------------

# Definieer de gewenste comparisons
selected_comparisons <- c("aCD3_aCD27_vs_aCD3", 
                          "aCD3_aCD28_vs_aCD3", 
                          "aCD3_a4_1BB_vs_aCD3", 
                          "aCD3_aCD27_aCD28_vs_aCD3", 
                          "aCD3_aCD28_a4_1BB_vs_aCD3")

# Filter de overlap matrix zodat alleen de geselecteerde comparisons overblijven (zowel rijen als kolommen)
overlap_matrix_filtered <- overlap_matrix[selected_comparisons, selected_comparisons]

# -------------------------------
# 1. Overlap-matrix herordenen
# -------------------------------
desired_order <- c("aCD3_aCD27_vs_aCD3", 
                   "aCD3_a4_1BB_vs_aCD3", 
                   "aCD3_aCD28_vs_aCD3",
                   "aCD3_aCD27_aCD28_vs_aCD3", 
                   "aCD3_aCD28_a4_1BB_vs_aCD3")

# Zorg dat de rijen en kolommen in dezelfde (handmatige) volgorde staan
overlap_matrix_filtered <- overlap_matrix_filtered[desired_order, desired_order]

# -------------------------------
# 2. Kleur- en schaalinstellingen
# -------------------------------
# Maak een vector van breaks van 0 t/m 100
my_breaks <- seq(0, 100, by = 1)

# Genereer een kleurenpalet (wit -> blauw)
my_colors <- colorRampPalette(c("white", "#FFC273"))(length(my_breaks) - 1)

# -------------------------------
# 3. Heatmap plotten (zonder clustering)
# -------------------------------
p <- pheatmap(
  overlap_matrix_filtered,
  color            = my_colors,
  breaks           = my_breaks,
  cluster_rows     = FALSE,     
  cluster_cols     = FALSE,
  display_numbers  = TRUE,      # Getallen in de cellen tonen
  number_format    = "%.1f",    # Afronden op 1 decimaal
  fontsize_number  = 12,        # Tekstgrootte voor de getallen in de cellen
  number_color     = "black",   # Kleur van de getallen in de cellen
  fontsize         = 12,        # Algemene tekstgrootte (legenda, titel)
  fontsize_row     = 12,        # Rijlabels
  fontsize_col     = 12,        # Kolomlabels
  angle_col        = 45,        # Draai kolomlabels 45 graden
  border_color     = "grey90",  # Lichte rand om de cellen
  main             = "Overlap (%) of DE genes between Comparisons compared to aCD3"
)

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/procent_heatmap_trans"
output_base <- paste0(output_dir, "procent_DE_transcriptomics")

# Sla op als PNG (hoge resolutie)
ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)


# Sla op als SVG
ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)


# Sla op als EPS (voor vectorformaat)
ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)


