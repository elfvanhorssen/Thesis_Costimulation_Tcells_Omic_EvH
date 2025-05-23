#venn CELL CYCLE 

library(biomaRt)
library(VennDiagram)
library(grid)

proteins_cd3_cd28_cd137 <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/2.unbiased_proteomics/results/top_genes/aCD3_aCD28_a4_1BB - aCD3_top_Genes_FC0.5.csv")
genes_cd3_cd28_cd137 <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/aCD3_vs_aCD3_aCD28_a4_1BB_DEGs.csv")


# make DEG and DEP sets 
deg_set <- unique(toupper(trimws(genes_cd3_cd28_cd137$geneName)))
dep_set <- unique(toupper(trimws(proteins_cd3_cd28_cd137$gene_id)))

# Use Ensembl for gene convertion
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# List some GO terms 
go_terms_cellcycle <- c("GO:0007049", # Cell cycle
                        "GO:0051301", # Cell division
                        "GO:0000278", # Mitotic cell cycle
                        "GO:0000086", # G2/M transition
                        "GO:0000075") # Cell cycle checkpoint

# Get protein code
cell_cycle_bm <- getBM(
  attributes = c("external_gene_name", "gene_biotype", "go_id"),
  filters = c("go", "biotype"),
  values = list(go_terms_cellcycle, "protein_coding"),
  mart = ensembl
)

# Get only unique genes 
cell_cycle_genes <- unique(toupper(trimws(cell_cycle_bm$external_gene_name)))

# Filter on the cell cycle gene set 
deg_cellcycle <- deg_set[deg_set %in% cell_cycle_genes]
dep_cellcycle <- dep_set[dep_set %in% cell_cycle_genes]

# Make Venn
venn.plot <- draw.pairwise.venn(
  area1 = length(deg_cellcycle),
  area2 = length(dep_cellcycle),
  cross.area = length(intersect(deg_cellcycle, dep_cellcycle)),
  category = c("Cell Cycle DEGs", "Cell Cycle DEPs"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 180),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE
)


output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Venn/"
output_base <- paste0(output_dir, "cell_cycle_cd3_cd28_cd137")

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

# export overlapping genes 
overlap_genes <- intersect(deg_coding, dep_coding)
write.csv(overlap_genes, paste0(output_base, "_overlap_cell_cycle_cd3_cd28_cd137.csv"), row.names = FALSE)