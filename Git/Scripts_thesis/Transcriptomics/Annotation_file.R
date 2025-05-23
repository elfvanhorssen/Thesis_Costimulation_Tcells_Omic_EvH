#Gene names for Ensembl ID

# Load library 
library(biomaRt)
library(dplyr)
library(tibble)

# open data 
annotation_file <- read.csv(gzfile("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/data/annotation.tsv.gz"), sep = "\t", header = TRUE)

#Select the right columns
annotation_file <- annotation_file %>% select(c(1:2))
