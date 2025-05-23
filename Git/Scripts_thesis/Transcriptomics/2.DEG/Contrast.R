#All DESeq2 comparisons in 1 list for further analysis

library(DESeq2)
groupscsp  <- c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")
results_listcsp_inverse  <- list()

#loop to fill results list

for (i in 1:(length(groupscsp )-1)){
  for (j in (i+1):length(groupscsp )){
    contrast <-c("condition", groupscsp [j], groupscsp [i])
    res_csp  <- results(dds_csp , contrast = contrast)
    results_listcsp_inverse [[paste(groupscsp [j], "vs", groupscsp [i], sep = "_")]] <- res_csp 
  }
}

#check if it went correct
names(results_listcsp_inverse)
summary(results_listcsp_inverse)

#save
write.csv(as.data.frame(results_listcsp_inverse), file = "I:\Research\TCR\2.Lab Members\Eralin van Horssen\co_stim_project\figures\deg\all/result_dds_list_21_unbiased_comparisons")
saveRDS(results_listcsp_inverse, file = "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/results_listcsp_inverse.rds")

