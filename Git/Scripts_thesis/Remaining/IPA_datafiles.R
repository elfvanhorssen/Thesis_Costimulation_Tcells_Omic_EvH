library(DESeq2)
library(openxlsx)  

# Define groups 
groupscsp <- c("Unstim", "aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")

# Make an empthy list 
results_listcsp <- list()

# Loop trough every comparison in the transcriptomic data 
for (i in 1:(length(groupscsp) - 1)) {
  for (j in (i + 1):length(groupscsp)) {
    
    # Define contrast for DESeq2
    contrast <- c("condition", groupscsp[i], groupscsp[j])
    res_csp <- results(dds_csp, contrast = contrast)
    
    # Give the comparison an name 
    comparison_name <- paste0(groupscsp[i], "_vs_", groupscsp[j])
    
    # make dataframe 
    res_csp_df <- as.data.frame(res_csp)
    
    # add gene names as a column
    res_csp_df$gene <- rownames(res_csp_df)
    
    # genes column first
    res_csp_df <- res_csp_df[, c("gene", setdiff(colnames(res_csp_df), "gene"))]
    
    # save in a list
    results_listcsp[[comparison_name]] <- res_csp_df
    
    # get files
    file_name <- paste0("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/result_dds_", comparison_name, ".xlsx")
    
    # write as excel 
    write.xlsx(res_csp_df, file = file_name)
  }
}
