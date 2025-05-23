library(DESeq2)
library(tibble)


ddds_vst <- rlog(dds_csp, blind = TRUE)
expr_matrix <- assay(ddds_vst)

#calculation of the z_score
z_scores <- t(scale(t(expr_matrix)))


#gene names 
z_scores <- as.data.frame(z_scores)
z_scores <- rownames_to_column(z_scores, var = "feature")
z_scores <- merge(z_scores, annotation_file, by = "feature")


z_scores$geneName <- make.unique(z_scores$geneName)
rownames(z_scores) <- z_scores$geneName

z_scores$geneName = NULL
z_scores$feature = NULL 

#good column names
metadata <- metadata
metadata$unique <- make.unique(as.character(metadata$condition))

# trim on spaces you cannot see
colnames(z_scores) <- trimws(colnames(z_scores))
rownames(metadata) <- trimws(rownames(metadata))

#Check if the column names of fragment_per_gene corresponds with the row names of metadata
all(colnames(z_scores) %in% rownames(metadata))

#same order design as fragments
metadata <- metadata[colnames(z_scores), ]

#check if they are in the same order
all(colnames(z_scores) == rownames(metadata))

new_colnames <- metadata$unique
names(new_colnames) <- rownames(metadata)
colnames(z_scores) <- new_colnames[colnames(z_scores)]


write.csv(z_scores, "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/figures/data/z_scores_trans.csv")