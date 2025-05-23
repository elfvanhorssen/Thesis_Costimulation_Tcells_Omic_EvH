library(dplyr)
library(tibble)
library(writexl)

#Function to get all the DEGs

get_all_degs <- function(res, annotation, comparison,
                         padj_cutoff = 0.01, l2fc_cutoff = 1) {
  # make dataframe, ensembl column and merge the annotation file
  df <- as.data.frame(res) %>%
    rownames_to_column("ensembl") %>%
    merge(annotation, by.x = "ensembl", by.y = "feature", all.x = TRUE)
  
  # columns need to be numeric and adjust NA values
  df$padj <- as.numeric(df$padj)
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  df$padj[is.na(df$padj)] <- 1
  
  # Filter on the parameters and select columns needed. 
    filter(padj < padj_cutoff,
           abs(log2FoldChange) > l2fc_cutoff) %>%
    select(ensembl,
           geneName,
           baseMean,
           log2FoldChange,
           lfcSE,
           stat,
           pvalue,
           padj)
  
  if (nrow(df_sig) == 0) {
    message("Geen DEGs voor ", comparison)
    return(NULL)
  }
  
  message("Aantal DEGs voor ", comparison, ": ", nrow(df_sig))
  return(df_sig)
}

# All the comparisons we need, to use in the loop. 
selected_comparisons <- c(
  "aCD3_aCD27_vs_aCD3",
  "aCD3_a4_1BB_vs_aCD3",
  "aCD3_aCD28_vs_aCD3",
  "aCD3_aCD27_aCD28_vs_aCD3",
  "aCD3_aCD28_a4_1BB_vs_aCD3"
)

degs_per_comparison <- list()
for (cmp in selected_comparisons) {
  df_degs <- get_all_degs(results_listcsp_inverse[[cmp]],
                          annotation_file,
                          cmp)
  if (!is.null(df_degs)) {
    degs_per_comparison[[ make.names(cmp) ]] <- df_degs
  }
}

# make an Excel, with a tab for every comparison
write_xlsx(degs_per_comparison, path = "all_DEGs_per_comparison.xlsx")
message("Klaar: alle DEGs weggeschreven naar all_DEGs_per_comparison.xlsx")

