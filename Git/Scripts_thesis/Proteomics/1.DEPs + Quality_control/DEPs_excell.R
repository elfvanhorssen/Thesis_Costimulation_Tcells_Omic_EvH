
# install.packages("writexl")

library(dplyr)
library(tibble)
library(writexl)

# function to get DE proteins
get_all_limma_deps <- function(res, comparison,
                               p_cutoff    = 0.05,
                               fc_cutoff   = 0.5) {
  

  df <- as.data.frame(res) %>%
    rownames_to_column("feature") %>%
    mutate(
      P.Value    = ifelse(is.na(P.Value), 1, P.Value)
    )
  
  # filter on parameters
  df_sig <- df %>%
    filter(
      P.Value < p_cutoff,
      abs(logFC) > fc_cutoff
    ) %>%
    select(
                
      protein, 
      gene_id,
      AveExpr,           
      logFC,             
      t,                 
      P.Value,           
      adj.P.Val          
    )
  
  if (nrow(df_sig) == 0) {
    message("Geen significante proteïnen voor ", comparison)
    return(NULL)
  }
  
  message("Aantal DE-proteïnen voor ", comparison, ": ", nrow(df_sig))
  return(df_sig)
}

# comparisons
selected_comparisons <- c(
  "aCD3_aCD28 - aCD3",
  "aCD3_aCD28_a4_1BB - aCD3",
  "aCD3_aCD27 - aCD3",
  "aCD3_aCD27_aCD28 - aCD3",
  "aCD3_a4_1BB - aCD3"
)

# list for degs
deps_per_comparison <- list()

for (cmp in selected_comparisons) {
  res_obj  <- results_list_log[[cmp]]
  df_deps  <- get_all_limma_deps(res_obj, cmp)
  if (!is.null(df_deps)) {
    sheet_name <- make.names(cmp)
    deps_per_comparison[[sheet_name]] <- df_deps
  }
}

# make excel with DEPs, every tab is an different comparison
output_file <- "all_DE_proteins_per_comparison.xlsx"
write_xlsx(deps_per_comparison, path = output_file)

message("Klaar! Excel-bestand weggeschreven naar:\n  ", output_file)

