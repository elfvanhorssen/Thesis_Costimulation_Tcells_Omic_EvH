# Laad benodigde packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

t_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_trans.csv")
rownames(t_z_score) <- t_z_score$X  
t_z_score$X <- NULL

p_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_proteomics.csv")
rownames(p_z_score) <- p_z_score$X  
p_z_score$X <- NULL

# Indien gewenst: bekijk even de eerste paar rijen
head(t_z_score)
head(p_z_score)

# Zorg dat je dezelfde genen/proteïnen overhoudt

# Haal de gewenste genen op uit df2 (stel dat de kolom "gene" heet)
selected_genes <- same_degs$same_genes
t_z_score_sub <- t_z_score[rownames(t_z_score) %in% selected_genes, ]
p_z_score_sub <- p_z_score[rownames(p_z_score) %in% selected_genes, ]


t_z_score_sub <- rownames_to_column(t_z_score_sub, var = "gene")
p_z_score_sub <- rownames_to_column(p_z_score_sub, var = "gene")


# -------------------------------
# 2. Data omzetten naar long format
# -------------------------------
# Hier veronderstellen we dat elke kolom een conditie is.
# We maken een lange tabel waarin we aangeven of het transcriptomics of proteomics betreft.

library(ggplot2)
library(dplyr)
library(tidyr)

# Voorbeeld: Stel dat in t_z_score_sub de conditienamen zo zijn:
# "unstimulated", "aCD3_aCD27_vs_aCD3", "aCD3_aCD28_vs_aCD3", etc.
# En in p_z_score_sub bijvoorbeeld:
# "UNSTIM", "CD3_CD27", "CD3_CD28", etc.
# Dan pas je ze eerst aan zodat ze overeenkomen.

# Transcriptomics data omzetten naar long format en conditienamen recoden:
trans_long <- t_z_score_sub %>%
  mutate(gene = rownames(t_z_score_sub)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "Condition",
    values_to = "z_score"
  ) %>%
  mutate(
    DataType = "Transcriptomics",
    # Pas hier de transcriptomics conditienamen aan zodat ze overeenkomen met de proteomics
    Condition = recode(Condition,
                       "unstimulated"         = "Unstim",
                       "aCD3_aCD27_vs_aCD3"   = "aCD3_aCD27",
                       "aCD3_aCD28_vs_aCD3"   = "aCD3_aCD28",
                       "aCD3_a4_1BB_vs_aCD3"   = "aCD3_a4_1BB",
                       "aCD3_aCD27_aCD28_vs_aCD3" = "aCD3_aCD27_aCD28",
                       "aCD3_aCD28_a4_1BB_vs_aCD3" = "aCD3_aCD28_a4_1BB"
    )
  )

# Proteomics data omzetten naar long format en conditienamen recoden:
prot_long <- p_z_score_sub %>%
  mutate(gene = rownames(p_z_score_sub)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "Condition",
    values_to = "z_score"
  ) %>%
  mutate(
    DataType = "Proteomics",
    # Pas hier de proteomics conditienamen aan zodat ze overeenkomen met de transcriptomics
    Condition = recode(Condition,
                       "UNSTIM"     = "Unstim",
                       "aCD3_aCD27" = "aCD3_aCD27",
                       "aCD3_aCD28" = "aCD3_aCD28",
                       "aCD3_a4_1BB" = "aCD3_a4_1BB",
                       "CD3_CD27_CD28" = "aCD3_aCD27_aCD28",
                       "CD3_CD28_4_1BB" = "aCD3_aCD28_a4_1BB"
    )
  )

# Nu beide datasets combineren op basis van 'gene' en 'Condition'
combined_scatter <- inner_join(
  trans_long,
  prot_long,
  by = c("gene", "Condition"),
  suffix = c("_trans", "_prot")
)

# Maak de scatter plot: x = proteomics z_score, y = transcriptomics z_score
scatter_plot <- ggplot(combined_scatter, aes(x = z_score_prot, y = z_score_trans)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Proteomics z-score",
    y = "Transcriptomics z-score",
    title = "Vergelijking Proteomics versus Transcriptomics Z-scores",
    subtitle = "Gecodeerde conditienamen zorgen voor een correcte merge"
  ) +
  theme_minimal()

print(scatter_plot)
