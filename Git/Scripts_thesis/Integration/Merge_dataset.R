#library 
library(dplyr)

#load z_scores data
t_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/figures/data/z_scores_trans.csv")
rownames(t_z_score) <- t_z_score$X  
t_z_score$X <- NULL

p_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/figures/data/z_scores_proteomics.csv")
rownames(p_z_score) <- p_z_score$X  
p_z_score$X <- NULL


#look at distribution 
t_z_score_matrix <- as.matrix(t_z_score)
p <- hist(t_z_score_matrix, breaks=50, main = "rna distribution z-score", col = "blue")
p

p_z_score_matrix <- as.matrix(p_z_score)
q <- hist(p_z_score_matrix, breaks=50, main = "proteomic distribution z-score", col = "blue")
q


boxplot(t_z_score, main = "Transcriptomics (Log2-transformed)", col = "lightblue")
boxplot(p_z_score, main = "Proteomics (Log2-transformed)", col = "lightgreen")


#Look at the intersect for the same amount of gene names 
t_rownames <- rownames(t_z_score)
p_rownames <- rownames(p_z_score)
both_df <- intersect(t_rownames, p_rownames)

#filter the genes for the heatmap 

t_filtered <- t_z_score[rownames(t_z_score) %in% both_df,]
p_filtered <- p_z_score[rownames(p_z_score) %in% both_df,]


#mergedata 

merged_df <- merge(t_filtered, p_filtered, by = "row.names")
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL

write.csv(merged_df, "merged.csv")


#example 

arobic_glycolysis <- c("ALDOA", "ENO1", "GAPDH","GPI", "HK1", "HK2", "LDHA", "PFKM", "PGAM1", "PKM", "PKMYT1", "SLC2A1")





