library(ggplot2)
library(dplyr)
library(ggrepel)
library(tibble)

#Function to create volcano plots 
plot_vulcano <- function(res, annotation, title, genes_highlight) {
  
  # 1. make from results a dataframe and make the rownames as a column calles ensemb
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "ensemb")
  
  # 2. Merge annotation data
  res_df <- res_df %>%
    left_join(annotation, by = c("ensemb" = "feature"))
  
  # 3. Change NA for 1
  res_df <- res_df %>%
    mutate(padj = if_else(is.na(padj), 1, padj))
  
  # 4. Catagorise for up and down. 
  res_df <- res_df %>%
    mutate(status = case_when(
      padj < 0.01 & log2FoldChange >  1 ~ "up",
      padj < 0.01 & log2FoldChange < -1 ~ "down",
      TRUE                              ~ "NS"
    )) %>%
    mutate(status = factor(status, levels = c("down", "NS", "up")))
  
  # 5. select genes that needs to be highlighted.  
  highlight <- res_df %>%
    filter(geneName %in% genes_highlight)
  
  # 6. Plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(color = "grey80", size = 2) +
    geom_point(data = highlight,
               aes(x = log2FoldChange, y = -log10(padj)),
               color = "red", size = 3) +
    geom_text_repel(data = highlight,
                    aes(label = geneName),
                    size = 3,
                    box.padding = 0.3,
                    point.padding = 0.3) +
    theme_classic(base_size = 14) +
    ggtitle(title) +
    labs(x = "log2FC", y = "-log10(p-adjusted)")
  
  return(p)
}

#genes of interest
genes_highlight <- c(
  "SELL",   "GLUT1",   "TNFRSF4", "TNFRSF9",
  "SLC7A1", "CDK1",    "MKI67",   "TOP2A",
  "CDKN1B", "MCM4",    "MCM6",    "RB1"
)

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Volcano/transcriptomics/same_axis/supp"

#loop to make the volcano plots 
for (comparison in names(results_listcsp_inverse)) {
  res_21 <- results_listcsp_inverse[[comparison]]
  
  plot <- plot_vulcano(
    res             = res_21,
    annotation      = annotation_file,
    title           = comparison,
    genes_highlight = genes_highlight
  )
  
  output_base <- file.path(output_dir, paste0(comparison, "_ANNO"))
  
  # save in different formats: png, svg and eps. 
  ggsave(paste0(output_base, ".png"), plot = plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_base, ".svg"), plot = plot, width = 8, height = 6)
  ggsave(paste0(output_base, ".eps"), plot = plot, width = 8, height = 6, device = cairo_ps)
}


