library(biomaRt)
library(DiffBind)
library(GenomicRanges)

# connect with ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# List of ensembl id's I want to investigate
ensembl_ids <- c("ENSG00000100453", "ENSG00000067225", "ENSG00000117394", "ENSG00000139514", "ENSG00000105281", "ENSG00000103257", "ENSG00000168003", "ENSG00000284967", "ENSG00000278540", "ENSG00000148773", "ENSG00000123358", "ENSG00000179388", "ENSG00000135476", "ENSG00000131747", "ENSG00000147536", "ENSG00000111276", "ENSG00000170312", "ENSG00000123374", "ENSG00000117399", "ENSG00000158402")  # Voeg je Ensembl IDs toe

# Data frame to save results 
results_df <- data.frame(
  ensembl_id = character(),
  hgnc_symbol = character(),
  contrast = integer(),
  overlap_count = integer(),
  p_value = numeric(),
  log2_fold_change = numeric(),
  stringsAsFactors = FALSE
)

# Loop trough ensembl ID's
for (ensembl_id in ensembl_ids) {
  # search for gene information
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "ensembl_gene_id",
    values = ensembl_id,
    mart = mart
  )
  
  if (nrow(gene_info) > 0) {
    # get gene information 
    hgnc_symbol <- gene_info$hgnc_symbol[1]  # take first value when there are more values available. 
    gene_gr <- GRanges(
      seqnames = paste0("chr", gene_info$chromosome_name[1]),
      ranges = IRanges(start = gene_info$start_position[1], end = gene_info$end_position[1])
    )
    
    # Loop trough contrast
    for (i in 1:length(diff$contrasts)) {
      diff_db <- dba.report(diff, contrast = i, bUsePval = TRUE, th = 0.5, fold = 1)
      overlaps <- subsetByOverlaps(diff_db, gene_gr)
      
      # see if genes overlap
      if (length(overlaps) > 0) {
        overlap_count <- length(overlaps)
        
        # get information from overlapping genes
        overlap_data <- as.data.frame(overlaps)
        overlap_pval <- overlap_data$p.value 
        overlap_log2fc <- overlap_data$Fold  
        
        # add overlap to result
        for (j in seq_len(nrow(overlap_data))) {
          results_df <- rbind(results_df, data.frame(
            ensembl_id = ensembl_id,
            hgnc_symbol = ifelse(hgnc_symbol == "", "NA", hgnc_symbol),  
            contrast = i,
            overlap_count = overlap_count,
            p_value = overlap_pval[j],
            log2_fold_change = overlap_log2fc[j],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  } else {
    cat(paste("Gene", ensembl_id, "not found in Ensembl/n"))
  }
}

# see results 
print(results_df)

# save
write.csv(results_df, "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/results/gene_results_with_overlaps2.csv", row.names = FALSE)
