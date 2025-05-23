#PCA plot, VSD: deleted noise, stabalised the variance without changing the data, n500 features are used. 
library(ggplot2)
library(DESeq2)
library(ggrepel)

vsdcsp <- vst(dds_csp, blind= FALSE)

#zet in de juiste volgorde 

vsdcsp$condition <- factor(ddscsp$condition, levels = c("Unstim", "aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD27_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28","aCD3_aCD28_a4_1BB","aCD3_aCD27_aCD28_a4_1BB"))

#use PCA and make the PCA plot
pcaData <- plotPCA(vsdcsp, intgroup = "condition", returnData = TRUE, ntop= 1000)
percentVar <- round(100* attr(pcaData, "percentVar"))

groupscsp  <- c("Unstim","aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")
subset_pcaData <- pcaData[pcaData$condition %in% groupscsp, ]
subset_pcaData$conditon <- factor(subset_pcaData, levels = groupscsp)

custom_colors <- c("Unstim" = "blue","aCD3" = "#000000", "aCD3_aCD27" = "#00A651", "aCD3_a4_1BB" = "orange", "aCD3_aCD28" = "#0072BC", "aCD3_aCD27_aCD28" = "#8A2BE2", "aCD3_aCD28_a4_1BB" = "#ED1C24")

pca_plot <- ggplot(subset_pcaData, aes(PC1, PC2, color=condition)) + 
  geom_point(size = 3) + 
  xlab(paste0( "PC1: ", percentVar[1], "%variance")) +
  ylab(paste0( "PC2: ", percentVar[2], "%variance")) +
  coord_fixed() +
  scale_color_manual(values = custom_colors) +
  theme_classic()

pca_plot

#save as png 
ggsave("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/1.PCA/PCA_plot_conditions.png", plot = pca_plot, dpi = 300)


library(ggplot2)

# Bereken aslimieten voor transcriptomics
x_min_t <- floor(min(subset_pcaData$PC1) / 10) * 10
x_max_t <- ceiling(max(subset_pcaData$PC1) / 10) * 10
y_min_t <- floor(min(subset_pcaData$PC2) / 10) * 10
y_max_t <- ceiling(max(subset_pcaData$PC2) / 10) * 10

# Definieer breaks en labels (om en om labelen, met 0 altijd zichtbaar)
x_breaks_t <- seq(x_min_t, x_max_t, by = 10)
y_breaks_t <- seq(y_min_t, y_max_t, by = 10)
x_labels_t <- ifelse(abs(x_breaks_t) %% 20 == 0, x_breaks_t, "")
y_labels_t <- ifelse(abs(y_breaks_t) %% 10 == 0, y_breaks_t, "")

# Maak de PCA plot met dezelfde lettertype-instelling als bij proteomics
pca_plot <- ggplot(subset_pcaData, aes(PC1, PC2, fill = condition)) + 
  geom_point(size = 4, shape = 21, color = "black", stroke = 1) +
  ggtitle("Global transcriptome") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial")  # Stel hier je gewenste lettertype in
  ) +
  scale_fill_manual(
    values = custom_colors,
    labels = c(
      "aCD3" = "aCD3",
      "aCD3_aCD27" = "aCD3 + aCD27",
      "aCD3_a4_1BB" = "aCD3 + aCD137",
      "aCD3_aCD28" = "aCD3 + aCD28",
      "aCD3_aCD27_aCD28" = "aCD3 + aCD27 + aCD28",
      "aCD3_aCD28_a4_1BB" = "aCD3 + aCD28 + aCD137"
    ),
    guide = guide_legend(override.aes = list(color = "black", shape = 21, stroke = 1.2))
  ) +
  guides(fill = guide_legend(title = NULL)) +
  scale_x_continuous(
    limits = c(x_min_t, x_max_t),
    breaks = x_breaks_t,
    labels = x_labels_t,
    minor_breaks = seq(x_min_t, x_max_t, by = 5)
  ) +
  scale_y_continuous(
    limits = c(y_min_t, y_max_t),
    breaks = y_breaks_t,
    labels = y_labels_t,
    minor_breaks = seq(y_min_t, y_max_t, by = 5)
  )

pca_plot


ggsave("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/PCA/PCA_plot_conditions.png", plot = pca_plot)







###################### diff? 

vsdcsp$condition <- factor(ddscsp$condition, levels = c("Unstim", "aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD27_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28","aCD3_aCD28_a4_1BB","aCD3_aCD27_aCD28_a4_1BB"))

#use PCA and make the PCA plot
pcaData <- plotPCA(vsdcsp, intgroup = "condition", returnData = TRUE, ntop= 2000)
percentVar <- round(100* attr(pcaData, "percentVar"))

groupscsp  <- c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")
subset_pcaData <- pcaData[pcaData$condition %in% groupscsp, ]
subset_pcaData$conditon <- factor(subset_pcaData, levels = groupscsp)

custom_colors <- c("aCD3" = "gray", "aCD3_aCD27" = "#00A651", "aCD3_a4_1BB" = "orange", "aCD3_aCD28" = "#0072BC", "aCD3_aCD27_aCD28" = "#8A2BE2", "aCD3_aCD28_a4_1BB" = "#ED1C24")

pca_plot <- ggplot(subset_pcaData, aes(PC1, PC2, fill=condition)) + 
  geom_point(size = 4, shape = 21, color = "black", stroke = 1) +
  ggtitle("Global transcriptome") +
  xlab(paste0( "PC1: ", percentVar[1], "%variance")) +
  ylab(paste0( "PC2: ", percentVar[2], "%variance")) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = custom_colors,
                    labels = c("aCD3" = "aCD3",
                               "aCD3_aCD27" = "aCD3 + aCD27",
                               "aCD3_a4_1BB" = "aCD3 + aCD137",
                               "aCD3_aCD28" = "aCD3 + aCD28",
                               "aCD3_aCD27_aCD28" = "aCD3 + aCD27 + aCD28",
                               "aCD3_aCD28_a4_1BB" = "aCD3 + aCD28 + aCD137"), guide = guide_legend(override.aes = list(color = "black", shape = 21, stroke = 1.2))) +
  guides(fill = guide_legend(title = NULL))



pca_plot

# Definieer de bestandsnamen en paden
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/PCA/"
output_base <- paste0(output_dir, "PCA_plot")

# Sla op als PNG (hoge resolutie)
ggsave(paste0(output_base, ".png"), plot = pca_plot, width = 8, height = 6, dpi = 300)

# Sla op als PDF
ggsave(paste0(output_base, ".pdf"), plot = pca_plot, width = 8, height = 6)

# Sla op als SVG
ggsave(paste0(output_base, ".svg"), plot = pca_plot, width = 8, height = 6)

# Sla op als TIFF (voor publicaties)
ggsave(paste0(output_base, ".tiff"), plot = pca_plot, width = 8, height = 6, dpi = 300, compression = "lzw")

# Sla op als EPS (voor vectorformaat)
ggsave(paste0(output_base, ".eps"), plot = pca_plot, width = 8, height = 6, device = cairo_ps)

ggsave("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/PCAplot_m.png", plot = pca_plot)

# geom_point(size = 4, shape = 21, color = "black", stroke = 1) +
# 
#   scale_x_continuous(
#     limits = c(x_min, x_max),
#     breaks = x_breaks,           # Hoofd-ticks elke 10 eenheden
#     labels = x_labels,           # Om en om een label, met 0 altijd zichtbaar
#     minor_breaks = seq(x_min, x_max, by = 5)
#   ) +
#   scale_y_continuous(
#     limits = c(y_min, y_max),
#     breaks = y_breaks,
#     labels = y_labels,
#     minor_breaks = seq(y_min, y_max, by = 5)
#   )

############################ 3D

# Installeer en laad benodigde pakketten
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("DESeq2", quietly = TRUE)) install.packages("DESeq2")

library(plotly)
library(ggplot2)
library(DESeq2)

# Zet condition factor levels
vsdcsp$condition <- factor(ddscsp$condition, levels = c("Unstim", "aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD27_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28","aCD3_aCD28_a4_1BB","aCD3_aCD27_aCD28_a4_1BB"))

# Haal de expressiematrix op en voer PCA uit
pca_matrix <- assay(vsdcsp)  # Haal de genormaliseerde expressiegegevens op
pca_result <- prcomp(t(pca_matrix), center = TRUE, scale. = TRUE)  # PCA uitvoeren

# Zet de PCA-resultaten om in een dataframe
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- colnames(vsdcsp)  # Sample-namen toevoegen
pca_df$Condition <- vsdcsp$condition  # Conditie-informatie toevoegen

# Bereken de variantiepercentages
percentVar <- round(100 * summary(pca_result)$importance[2, 1:3], 1)

# Definieer de condities die je wilt plotten
groupscsp  <- c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")
subset_pca_df <- pca_df[pca_df$Condition %in% groupscsp, ]

# Definieer aangepaste kleuren
custom_colors <- c("aCD3" = "#000000", "aCD3_aCD27" = "#00A651", "aCD3_a4_1BB" = "orange", "aCD3_aCD28" = "#0072BC", "aCD3_aCD27_aCD28" = "#8A2BE2", "aCD3_aCD28_a4_1BB" = "#ED1C24")

# Maak een interactieve 3D PCA-plot met plotly
fig <- plot_ly(subset_pca_df,
               x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~Condition, colors = custom_colors,
               text = ~paste("Sample:", Sample, "<br>Condition:", Condition),  # Hovertekst
               type = "scatter3d", mode = "markers", marker = list(size = 6, opacity = 0.8))

# Voeg labels en titel toe
fig <- fig %>% layout(title = "3D PCA: Global Transcriptome",
                      scene = list(xaxis = list(title = paste0("PC1: ", percentVar[1], "% variance")),
                                   yaxis = list(title = paste0("PC2: ", percentVar[2], "% variance")),
                                   zaxis = list(title = paste0("PC3: ", percentVar[3], "% variance"))))

# Toon plot
fig
