#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DiffBind")

 library(DiffBind)
# 
# #set working directory
 setwd("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac2/atac2/data")
# 
# # read in sampleSheet
# diff <- dba(sampleSheet="I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac2/atac2/data/sample_input_eralin2.csv", minOverlap=2, scoreCol=5)
# diff
# 
# # count reads in peaks
# diff <- dba.count(diff, summits=250)
# diff
# 
# saveRDS(diff, file = "dba_counts2.rds")
 

 
# # set up contrast
# diff <- dba.contrast(diff, categories=DBA_TREATMENT, minMembers=2)
# diff

diff <- readRDS("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac2/data/dba_counts2.rds")
diff 

# Make every contrast possible. 
# Get unique treatment 
treatment <- unique(diff$class["Treatment",])

#get all possible combinations 
contrast_combination <- combn(treatment, 2)

# loop trough combinations and make contrasts. 
for (i in 1:ncol(contrast_combination)) { 
  group1 <- contrast_combination[,i][1]
  group2 <- contrast_combination[,i][2]
  
  diff <- dba.contrast(
    diff, 
    group1 = diff$masks[[group1]],
   group2 = diff$masks[[group2]],
   name1 = group1,
   name2 = group2
 )
  
}

dba.show(diff, bContrasts = TRUE)
# perform edgeR- and DESeq2-based differential analysis
diff <- dba.analyze(diff, method=DBA_ALL_METHODS)
diff

# diff is analyzed 
# diff <- dba.analyze(diff, method=DBA_ALL_METHODS)

# Show number of contrast 
dba.show(diff, bContrasts = TRUE)

# Loop trough contrasts 
for (contrast_index in 1:length(diff$contrasts)) {
  # get the mask
  group1 <- diff$contrasts[[contrast_index]]$name1
  group2 <- diff$contrasts[[contrast_index]]$name2
  
  # Combine the nams of the groups 
  base_filename <- paste0(group1, "_vs_", group2)
  
  # get EdgeR-results
  edger_results <- dba.report(diff, contrast=contrast_index, method=DBA_EDGER)
  if (!is.null(edger_results)) {
    edger_filename <- paste0("edger_results_", base_filename, ".csv")
    write.csv(edger_results, file=edger_filename, row.names=FALSE)
    cat("EdgeR results saved to:", edger_filename, "/n")
  } else {
    cat("No EdgeR results for contrast:", base_filename, "/n")
  }
  
  # get DESeq2-results
  deseq_results <- dba.report(diff, contrast=contrast_index, method=DBA_DESEQ2)
  if (!is.null(deseq_results)) {
    deseq_filename <- paste0("deseq_results_", base_filename, ".csv")
    write.csv(deseq_results, file=deseq_filename, row.names=FALSE)
    cat("DESeq2 results saved to:", deseq_filename, "/n")
  } else {
    cat("No DESeq2 results for contrast:", base_filename, "/n")
  }
}



# save 
save.image(paste("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac2/atac2/data","diffbind_atac_h301_aml2.RData", sep=""))

# write csv file of differentially bound regions
diff.DB <- dba.report(diff, bUsePval = TRUE, th=0.5, fold = 1)
write.csv(diff.DB, "diffbind_atac_h301_aml2.csv")

dba.plotPCA(diff)


