#making the full DEGs list start with an empthy list 
library(dplyr)
library(DESeq2)
library(tidyverse)


top_genes_list_full <- list()

select_top_genes <- function(res, annotation, comparison) {
  
  #replace NAs for a high value so this genes will be placed at the bottom of the plot. 
  res$padj[is.na(res$padj)] <- 1
  
  res <- as.data.frame(res)
  res <- res %>% rownames_to_column(var = "ensemb")
  res <- merge(res, annotation, by.x = "ensemb", by.y = "feature", all.x = TRUE)
 
  
  #results in df and select top genes with filter on adj and l2fc
  filtered_res <- res[res$padj <0.01 & abs(res$log2FoldChange) >1, ]
  
  #controle if the comparison has significant genes 
  if (nrow(filtered_res) == 0) {
    cat("no genes found for", comparison, "\n")
    return(NULL)
  } else {
    cat("top genes for", comparison, ":\n")
    print(rownames(filtered_res))
    return(data.frame(Comparison = comparison, Gene = rownames(filtered_res), geneName = filtered_res$geneName))
    
  }
}

#loop for every comparison
print(names(results_listcsp))

selected_comparisons <- c("aCD3_vs_aCD3_aCD27", "aCD3_vs_aCD3_a4_1BB", "aCD3_vs_aCD3_aCD28", "aCD3_vs_aCD3_aCD27_aCD28", "aCD3_vs_aCD3_aCD28_a4_1BB")

for (comparison in selected_comparisons) {
  res_top_genes <- results_listcsp[[comparison]]
  top_genes <- select_top_genes(res_top_genes, annotation_file, comparison)
  if(!is.null(top_genes)) {
    top_genes_list_full[[comparison]] <- top_genes
  }
}

#save and check 
if (length(top_genes_list_full)>0){
  full_list <- do.call(rbind, top_genes_list_full)
  write.csv(full_list, file = "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/full_list_top_genes_costim.csv", row.names = FALSE)
} else {
  cat("the list is empthy")
}


#make a list were the genes are unique and save
unique_full_list <- full_list[!duplicated(full_list$Gene), ]
write.csv(unique_full_list, file = "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/unique_full_list_top_genes_costim.csv", row.names = FALSE )