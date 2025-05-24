#packages
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")

library(biomaRt)
library(VennDiagram)
library(grid)

# load data
proteins <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/2.unbiased_proteomics/results/top_genes/unique_full_list_top_GenesFC05_costim.csv")
genes <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/unique_full_list_top_genes_costim.csv")

# make unique
deg_set <- unique(toupper(trimws(genes$geneName)))
dep_set <- unique(toupper(trimws(proteins$Gene)))

# get the protein coding genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

coding_genes_df <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

#make unique list of coding genes
coding_genes <- unique(toupper(trimws(coding_genes_df$external_gene_name)))

# Filter on coding genes
deg_coding <- deg_set[deg_set %in% coding_genes]
dep_coding <- dep_set[dep_set %in% coding_genes]

# Make Venn
venn.plot <- draw.pairwise.venn(
  area1 = length(deg_coding),
  area2 = length(dep_coding),
  cross.area = length(intersect(deg_coding, dep_coding)),
  category = c("Coding DEGs", "Coding DEPs"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 180),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE
)


output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Venn/"
output_base <- paste0(output_dir, "new_coding_degs_vs_deps")

#  save
svg(paste0(output_base, ".svg"), width = 6, height = 6)
grid.draw(venn.plot)
dev.off()


cairo_ps(paste0(output_base, ".eps"), width = 6, height = 6, fallback_resolution = 300)
grid.draw(venn.plot)
dev.off()

png(paste0(output_base, ".png"), width = 1200, height = 1000, res = 300)
grid.draw(venn.plot)
dev.off()

# export the overlapping genes
overlap_genes <- intersect(deg_coding, dep_coding)
write.csv(overlap_genes, paste0(output_base, "_overlap_genes.csv"), row.names = FALSE)
