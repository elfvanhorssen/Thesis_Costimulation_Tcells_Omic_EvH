### 1️⃣ Vereiste packages laden
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)

# Load data
deseq_res <- readRDS("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/results_listcsp_inverse.rds")



### 3️⃣ Up/Downregulated genen tellen per comparison
deg_counts <- data.frame(Comparison = character(), Up = integer(), Down = integer(), stringsAsFactors = FALSE)

# Loop door elke vergelijking in de lijst
for (comparison in names(deseq_res)) {
  res <- deseq_res[[comparison]]
  
  # Genen categoriseren als Up, Down of Unchanged
  up_count <- sum(res$padj < 0.01 & res$log2FoldChange > 1, na.rm = TRUE)
  down_count <- sum(res$padj < 0.01 & res$log2FoldChange < -1, na.rm = TRUE)
  
  # Opslaan in dataframe
  deg_counts <- rbind(deg_counts, data.frame(Comparison = comparison, Up = up_count, Down = down_count))
}



library(tidyr)
library(ggplot2)

# Stel, je dataframe 'deg_counts' bevat de kolommen: Comparison, Up en Down

# Zet de data om van wide naar long formaat:
deg_long <- pivot_longer(deg_counts,
                         cols = c("Up", "Down"),
                         names_to = "Regulation",
                         values_to = "Count")


selected_comparisons <- c("aCD3_aCD27_vs_aCD3", "aCD3_aCD28_vs_aCD3", "aCD3_a4_1BB_vs_aCD3", "aCD3_aCD27_aCD28_vs_aCD3", "aCD3_aCD28_a4_1BB_vs_aCD3")

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter alleen de gewenste comparisons
selected_comparisons <- c("aCD3_aCD27_vs_aCD3", 
                          "aCD3_aCD28_vs_aCD3", 
                          "aCD3_a4_1BB_vs_aCD3", 
                          "aCD3_aCD27_aCD28_vs_aCD3", 
                          "aCD3_aCD28_a4_1BB_vs_aCD3")

df_filtered <- deg_long %>% 
  filter(Comparison %in% selected_comparisons)

# Stel de gewenste volgorde in voor de Comparison variabele
desired_order <- c("aCD3_aCD27_vs_aCD3", 
                   "aCD3_a4_1BB_vs_aCD3", 
                   "aCD3_aCD28_vs_aCD3",
                   "aCD3_aCD27_aCD28_vs_aCD3", 
                   "aCD3_aCD28_a4_1BB_vs_aCD3")

df_filtered$Comparison <- factor(df_filtered$Comparison, levels = desired_order)

# Indien gewenst, stel ook de volgorde van de Regulation in (bijvoorbeeld eerst Up, dan Down)
df_filtered$Regulation <- factor(df_filtered$Regulation, levels = c("Up", "Down"))

# Maak de stacked bar chart met de gefilterde data
plot1 <- ggplot(df_filtered, aes(x = Count, y = Comparison, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Up" = "salmon", "Down" = "lightblue")) +
  labs(title = "Differential expressed Up- and Downregulated genes per comparison, padj < 0.01. logfc abs 1",
       x = "genes", y = "Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot1)


library(dplyr)
library(ggplot2)

# 1. Maak 'Down' negatief
df_diverging <- df_filtered %>%
  mutate(Count = ifelse(Regulation == "Down", -Count, Count))

# 2. Plotten als 'diverging bar chart'
p <- ggplot(df_diverging, aes(x = Count, y = Comparison, fill = Regulation)) +
  #geom_col() +
  geom_col(width = 0.5) + 
  # Verticale lijn op x=0 zodat duidelijk is waar de scheiding tussen Up/Down zit
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = c("Up" = "salmon", "Down" = "lightblue")) +
  labs(title = "DE genes padj < 0.01 & res$log2FoldChange > 1 compared to aCD3",
       x = "Number of genes",
       y = "Comparison") +
  theme_classic()

p

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/"
output_base <- paste0(output_dir, "DE_transcriptomics")

# Sla op als PNG (hoge resolutie)
ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)


# Sla op als SVG
ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)


# Sla op als EPS (voor vectorformaat)
ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)

