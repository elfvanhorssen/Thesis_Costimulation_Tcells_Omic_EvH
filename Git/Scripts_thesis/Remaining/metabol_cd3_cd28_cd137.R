library(biomaRt)
library(VennDiagram)
library(grid)

proteins_cd3_cd28_cd137 <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/2.unbiased_proteomics/results/top_genes/aCD3_aCD28_a4_1BB - aCD3_top_Genes_FC0.5.csv")
genes_cd3_cd28_cd137 <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/aCD3_vs_aCD3_aCD28_a4_1BB_DEGs.csv")


# 📌 Je eigen DEGs en DEPs (schoongemaakt)
deg_set <- unique(toupper(trimws(genes_cd3_cd28_cd137$geneName)))
dep_set <- unique(toupper(trimws(proteins_cd3_cd28_cd137$gene_id)))

# 📡 Connectie met Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


# 🧬 GO-terms voor metabolism (algemeen + subsets)
go_terms_metabolism <- c("GO:0008152",  # Metabolic process
                         "GO:0006091",  # Generation of precursor metabolites and energy
                         "GO:0006096",  # Glycolytic process
                         "GO:0006629",  # Lipid metabolic process
                         "GO:0005975",  # Carbohydrate metabolic process
                         "GO:0006082",
                         "GO:0045730",
                         "GO:1901160",
                         "GO:0090345",
                         "GO:0019748",
                         "GO:0001887",
                         "GO:1902421",
                         "GO:0009812",
                         "GO:0018942",
                         "GO:0018930",
                         "GO:1901275",
                         "GO:0070085",
                         "GO:1990845",
                         "GO:0140651",
                         "GO:0019222",
                         "GO:1901440",
                         "GO:0050898",
                         "GO:0042620",
                         "GO:0042430",
                         "GO:0010191",
                         "GO:0006518",
                         "GO:0008218",
                         "GO:0043603",
                         "GO:0061984",
                         "GO:0072593",
                         "GO:0032963",
                         "GO:0009308"
                         )  # Organic acid metabolic process

# 🔍 Ophalen van protein-coding genen met deze GO-termen
metabolism_bm <- getBM(
  attributes = c("external_gene_name", "gene_biotype", "go_id"),
  filters = c("go", "biotype"),
  values = list(go_terms_metabolism, "protein_coding"),
  mart = ensembl
)

# 🧼 Unieke genen opschonen
metabolism_genes <- unique(toupper(trimws(metabolism_bm$external_gene_name)))

# 🔎 Filter jouw sets
deg_metabolism <- deg_set[deg_set %in% metabolism_genes]
dep_metabolism <- dep_set[dep_set %in% metabolism_genes]

# 📊 Venn plot
venn.plot <- draw.pairwise.venn(
  area1 = length(deg_metabolism),
  area2 = length(dep_metabolism),
  cross.area = length(intersect(deg_metabolism, dep_metabolism)),
  category = c("Metabolism DEGs", "Metabolism DEPs"),
  fill = c("lightgreen", "orange"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 180),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE
)

# Outputpad (pas eventueel aan)
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Venn/"
output_base <- paste0(output_dir, "metabol_cd3_cd28_cd137")

#  Opslaan als SVG
svg(paste0(output_base, ".svg"), width = 6, height = 6)
grid.draw(venn.plot)
dev.off()

# Opslaan als EPS (vector)
cairo_ps(paste0(output_base, ".eps"), width = 6, height = 6, fallback_resolution = 300)
grid.draw(venn.plot)
dev.off()

# Opslaan als PNG (raster)
png(paste0(output_base, ".png"), width = 1200, height = 1000, res = 300)
grid.draw(venn.plot)
dev.off()

# exporteer overlappende coding genen als CSV
overlap_genes <- intersect(deg_coding, dep_coding)
write.csv(overlap_genes, paste0(output_base, "_overlap_metabol_cd3_cd28_cd137.csv"), row.names = FALSE)
