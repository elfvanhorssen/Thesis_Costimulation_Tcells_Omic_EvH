library(ggplot2)


#open tz file with csv 
fragments_per_gene_check <- read.csv(gzfile("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/data/all_samples.fragments_per_gene.gz"), sep = "\t", header = TRUE)

# make the first column of the fragments_per_gene as row names
row.names(fragments_per_gene_check) <- fragments_per_gene_check[, 1]
#fragments_per_gene_check <- fragments_per_gene_check[, -1]


#calculate library size 
library_size <- colSums(fragments_per_gene_check)
library_size


ggplot(data.frame(sample = names(library_size), LibrarySize = library_size), aes(x=colnames(fragments_per_gene_check), y=LibrarySize)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Library Size Distribution") +
  ylab("Total Counts") +
  xlab("Samples")+
  theme(axis.text.x = element_text(angle = 45, hjust =1))


#check zero counts 
zero_counts_per_sample <- colSums(fragments_per_gene_check == 0)

#take percentages 
percentage_zero_counts_per_sample <- zero_counts_per_sample/ nrow(fragments_per_gene_check) * 100 

#plot
ggplot(data.frame(sample = names(percentage_zero_counts_per_sample), ZeroCounts = percentage_zero_counts_per_sample), aes(x=colnames(fragments_per_gene_check), y=ZeroCounts)) + 
  geom_bar(stat = "identity", fill="steelblue") +
  theme_minimal() +
  ggtitle("Percentage of zero counts per sample") +
  ylab("Percentage of zero counts") +
  xlab("Samples")+
  theme(axis.text.x = element_text(angle = 45, hjust =1))


#distirbution of gene counts per sample
distribution_gene_counts_per_sample <- ggplot(stack(as.data.frame(fragments_per_gene_check)), aes(x=ind, y=values)) +
  geom_boxplot()+
  theme_minimal() +
  ggtitle("Gene count distribution across samples") +
  ylab("Counts per gene") +
  xlab("Samples")+
  theme(axis.text.x = element_text(angle = 45, hjust =1))



