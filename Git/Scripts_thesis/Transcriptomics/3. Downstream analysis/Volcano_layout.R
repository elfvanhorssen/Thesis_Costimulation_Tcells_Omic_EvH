library(ggplot2)
library(dplyr)
library(ggrepel)
library(tibble)

#function to make volcano plots 
plot_vulcano <- function(res, annotation, title) {
  # results need to be an dataframe and add ensembl column. 
  res <- as.data.frame(res)
  res <- res %>% rownames_to_column(var = "ensemb")
  res <- merge(res, annotation, by.x = "ensemb", by.y = "feature", all.x = TRUE)
  res$padj[is.na(res$padj)] <- 1
  
  res_df <- as.data.frame(res)
  
  # Make column with status: upregulated, downregulated and not significant (NS)
  res_df$status <- "NS"
  res_df$status[res_df$padj < 0.01 & res_df$log2FoldChange > 1] <- "up"
  res_df$status[res_df$padj < 0.01 & res_df$log2FoldChange < -1] <- "down"
  res_df$status <- factor(res_df$status, levels = c("down", "NS", "up"))
  
  # Selected top 8 genes for labeling
  filtered_res <- res_df[res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1, ]
  top_genes <- head(filtered_res[order(filtered_res$padj), ], 8)
  
  # Count up and downregulated genes. 
  count_down <- sum(res_df$status == "down", na.rm = TRUE)
  count_up   <- sum(res_df$status == "up", na.rm = TRUE)
  
  # make titels for the plots 
  cond_label <- sub("_vs_.*", "", title)
  
  # make the de volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(fill = status), shape = 21, size = 2.5, color = "black", stroke = 0.3) +
    theme_classic(base_size = 14) +
    ggtitle(title) +
    scale_fill_manual(values = c("down" = "lightblue", "up" = "salmon")) +
    labs(x = "l2fc", y = "-log10 adjusted p-value") +
    
    # make the axis for every plot the same 
    scale_x_continuous(limits = c(-15, 15), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 80), expand = expansion(0)) +
    
    # help lines 
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "gray") +
    
    # Label the top 8 genes 
    geom_text_repel(
      data = top_genes,
      aes(label = geneName),
      size = 2.5,
      box.padding = 0.4,
      point.padding = 0.5,
      max.overlaps = 20
    ) +
    
    # add text for what is shown in the data 
    annotate(
      "text", x = -Inf, y = Inf,
      label = paste(cond_label, "Down"),
      hjust = -0.1, vjust = 1.2,  
      size = 4, fontface = "bold"
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste(cond_label, "Up"),
      hjust = 1.1, vjust = 1.2,
      size = 4, fontface = "bold"
    ) +
    
    # add text with the up and down-regulated counts
    annotate(
      "text", x = -Inf, y = -Inf,
      label = paste("n =", count_down),
      hjust = -0.1, vjust = -0.2,
      size = 4, color = "lightblue"
    ) +
    annotate(
      "text", x = Inf, y = -Inf,
      label = paste("n =", count_up),
      hjust = 1.1, vjust = -0.2,
      size = 4, color = "salmon"
    ) +
    
    # labels needs to be in the figure 
    coord_cartesian(clip = "off")
  
  return(p)
}

# loop
for (comparison in names(results_listcsp_inverse)) {
  res_21 <- results_listcsp_inverse[[comparison]]
  plot <- plot_vulcano(res_21, annotation_file, comparison)
  output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Volcano/transcriptomics/same_axis/supp/"
  output_base <- paste0(output_dir, comparison, "vol_same_axis_trans")
  
  # Saving of the figures
  ggsave(paste0(output_base, ".png"), plot = plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_base, ".svg"), plot = plot, width = 8, height = 6)
  ggsave(paste0(output_base, ".eps"), plot = plot, width = 8, height = 6, device = cairo_ps)
  
  
}
