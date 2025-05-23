#calculate z_scores 
library(tibble)
library(preprocessCore)
annotation <- anno


z_score <- protein_data_logged
z_score <- as.data.frame(z_score)
z_score <- rownames_to_column(z_score, var = "protein")
z_score <- merge(z_score, annotation, by = "protein" )



z_score$gene_id[is.na(z_score$gene_id)] <- "unkown"
z_score$gene_id <- make.unique(z_score$gene_id)
sum(is.na(z_score$gene_id))
rownames(z_score) <- z_score$gene_id
z_score$protein <- NULL
z_score$gene_id <- NULL
z_score <- t(scale(t(z_score)))


# z_score <- protein_log_gene_done
# z_score$Protein <- NULL
# z_score$Gene.Symbol <- NULL
# z_score <- 
# z_score <- t(scale(t(z_score)))


write.csv(z_score, "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/figures/data/z_scores_proteomics.csv")
