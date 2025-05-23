# install.packages("tidyverse")
# library(BiocManager)
# BiocManager:: install("airway")
# BiocManager:: install("DESeq2")
#BiocManager:: install("apeglm")
# library(tidyverse)
# library(airway)
# library(dplyr)
# library(apeglm)
library(pheatmap)
library(DESeq2)
library(ggrepel)

#open tz file with csv 
fragments_per_gene <- read.csv(gzfile("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/data/all_samples.fragments_per_gene.gz"), sep = "\t", header = TRUE)

# make the first column of the fragments_per_gene as row names
row.names(fragments_per_gene) <- fragments_per_gene[, 1]
fragments_per_gene <- fragments_per_gene[!grepl("__alignment_not_unique|__ambiguous|__no_feature|__not_aligned|__too_low_aQual", fragments_per_gene$feature), ]
#fragments_per_gene <- fragments_per_gene[, -1]

# open design matrix 
design_matrix <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/data/design_matrix_raw_collection_day.tsv", sep = "\t", header = TRUE)

# make the first column of the design_matrix as row names
row.names(design_matrix) <- design_matrix[, 1]
#design_matrix <- design_matrix[, -1]

# trim on spaces you cannot see
colnames(fragments_per_gene) <- trimws(colnames(fragments_per_gene))
rownames(design_matrix) <- trimws(rownames(design_matrix))

#Check if the column names of fragment_per_gene corresponds with the row names of design_matrix
all(colnames(fragments_per_gene) %in% rownames(design_matrix))

#same order design as fragments
design_matrix <- design_matrix[colnames(fragments_per_gene), ]

#check if they are in the same order
all(colnames(fragments_per_gene) == rownames(design_matrix))

#make design from character a vector
design_matrix$condition <- as.factor(design_matrix$condition)

#start deseq
ddscsp <- DESeqDataSetFromMatrix(countData = fragments_per_gene,
                              colData = design_matrix,
                              design = ~ condition + sex)

ddscsp

# pre-filtering, remove rows with low gene counts, keep rows that have at least 10 reads. 
keep <- rowSums(counts(ddscsp)) >= 10 
ddscsp <- ddscsp[keep,]
print(ddscsp)

#set the factor level, compare unstim with others
ddscsp$condition <- relevel(ddscsp$condition, ref = "Unstim")
ddscsp

#run deseq 
dds_csp <- DESeq(ddscsp)

#save results 
res_csp <- results(dds_csp, alpha = 0.01)
res_csp
