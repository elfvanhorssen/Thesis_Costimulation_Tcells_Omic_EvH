library(openxlsx)
library(VennDiagram)

# what do we want to show? comparisons and there overlapping degs (up and down). 

#load the data 


unstim_vs_cd3_proteins <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/2.unbiased_proteomics/results/top_genes/Unstim - aCD3_top_Genes.csv")

unstim_vs_cd3_genes <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/Unstim_vs_aCD3_DEGs.csv")

unstim_vs_cd3_deps <- unstim_vs_cd3_proteins$gene_id
unstim_vs_cd3_degs <- unstim_vs_cd3_genes$geneName

#zelfde naamgeving in de lijten (met hoofdletters en trimm etcc)

unstim_vs_cd3_deps <- toupper(trimws(unstim_vs_cd3_deps))
unstim_vs_cd3_degs <- toupper(trimws(unstim_vs_cd3_degs))


#make sets 

dep_set <- unique(unstim_vs_cd3_deps[!is.na(unstim_vs_cd3_deps)])
deg_set <- unique(unstim_vs_cd3_degs[!is.na(unstim_vs_cd3_degs)])

#make venn 

plot <- venn.diagram(
  x = list(DEGs = deg_set, DEPs = dep_set),
  category.names = c("DEGs", "DEPs"),
  filename = ("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/3.unbiased_indirect_integration/results/venn2/Venn unstim_vs_cd3.png"),
  output = TRUE,
  fill = c("lightblue", "lightpink"),
  alpha = 0.6,
  cex = 1,
  height = 5000,
  width = 5000,
  margin = 0.1,
  scaled = FALSE
)

venn.plot <- draw.pairwise.venn(
  area1 = length(deg_set),
  area2 = length(dep_set),
  cross.area = length(intersect(deg_set, dep_set)),
  catagory = c("deg_set", "dep_set"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  scaled = TRUE
)



png("Venn unstim_vs_cd3.png")
grid.draw(venn.plot)
dev.off()
