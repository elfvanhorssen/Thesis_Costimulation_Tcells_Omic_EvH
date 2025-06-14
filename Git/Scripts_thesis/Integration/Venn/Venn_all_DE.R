library(openxlsx)
library(VennDiagram)
library(grid)

#load the data 
proteins <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/2.unbiased_proteomics/results/top_genes/unique_full_list_top_GenesFC05_costim.csv")
genes <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/4.top_genes/unique_full_list_top_genes_costim.csv")

#unique deg and deps set 
deg_set <- unique(toupper(trimws(genes$geneName)))
dep_set <- unique(toupper(trimws(proteins$Gene)))

# Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(deg_set),
  area2 = length(dep_set),
  cross.area = length(intersect(deg_set, dep_set)),
  category = c("Transcriptomics (DEG)", "Proteomics (DEP)"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 2,            
  cat.cex = 2,        
  cat.pos = c(0, 180),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE       
)

# output directiory 
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/Venn/"
output_base <- paste0(output_dir, "all_degs_vs_deps")

# save
svg(paste0(output_base, ".svg"), width = 6, height = 6)
grid.draw(venn.plot)
dev.off()

cairo_ps(filename = paste0(output_base, ".eps"), width = 6, height = 6, fallback_resolution = 300)
grid.draw(venn.plot)
dev.off()

png(paste0(output_base, ".png"), width = 1200, height = 1000, res = 300)
grid.draw(venn.plot)
dev.off()
