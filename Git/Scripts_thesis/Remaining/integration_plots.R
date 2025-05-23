#code for density 2D plot 

library(dplyr)
library(tidyr)
library(ggplot2)

t_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_trans.csv")
p_z_score <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/R/data/z_scores_proteomics.csv")


# Stel je hebt een dataframe 'same_genes' met een kolom 'Gene' waarin de genen staan.
#genes_of_interest <- same_degs$same_genes


# Schoonmaken en uniek maken (toupper + trimws)


deg_coding_venn <- deg_set[deg_set %in% coding_genes]
dep_coding_venn <- dep_set[dep_set %in% coding_genes]

deg_coding_venn_use <- unique(toupper(trimws(deg_coding_venn)))
dep_coding_venn_use<- unique(toupper(trimws(dep_coding_venn)))
#look at which are the same

genes_of_interest <- intersect(deg_coding_venn_use, dep_coding_venn_use)
length(genes_of_interest)


# 1. Transcriptomics naar long + filter op genen
t_long <- t_z_score %>%
  filter(X %in% genes_of_interest) %>%
  pivot_longer(
    cols = -X,  # Alle kolommen behalve de gennaam
    names_to = "Condition_t",
    values_to = "Transcriptomics_z"
  ) 


library(dplyr)

t_long <- t_long %>%
  mutate(condition_cleaned = recode(Condition_t,
                                    "Unstim" = "unstim",
                                    "Unstim.1" = "unstim",
                                    "Unstim.2" = "unstim",
                                    "Unstim.3" = "unstim",
                                    "Unstim.4" = "unstim",
                                    "aCD3.1"   = "aCD3",
                                    "aCD3.2"   = "aCD3",
                                    "aCD3.3"   = "aCD3",
                                    "aCD3.4"   = "aCD3",
                                    "aCD3_aCD27.1"   = "aCD3_aCD27",
                                    "aCD3_aCD27.2"   = "aCD3_aCD27",
                                    "aCD3_aCD27.3"   = "aCD3_aCD27",
                                    "aCD3_aCD27.4"   = "aCD3_aCD27",
                                    "aCD3_a4_1BB.1"   = "aCD3_a4_1BB",
                                    "aCD3_a4_1BB.2"   = "aCD3_a4_1BB",
                                    "aCD3_a4_1BB.3"   = "aCD3_a4_1BB",
                                    "aCD3_a4_1BB.4"   = "aCD3_a4_1BB",
                                    "aCD3_aCD28.1"   = "aCD3_aCD28",
                                    "aCD3_aCD28.2"   = "aCD3_aCD28",
                                    "aCD3_aCD28.3"   = "aCD3_aCD28",
                                    "aCD3_aCD28.4"   = "aCD3_aCD28",
                                    "aCD3_aCD27_aCD28.1"   = "aCD3_aCD27_aCD28",
                                    "aCD3_aCD27_aCD28.2"   = "aCD3_aCD27_aCD28",
                                    "aCD3_aCD27_aCD28.3"   = "aCD3_aCD27_aCD28",
                                    "aCD3_aCD27_aCD28.4"   = "aCD3_aCD27_aCD28",
                                    "aCD3_aCD28_a4_1BB.1"   = "aCD3_aCD28_a4_1BB",
                                    "aCD3_aCD28_a4_1BB.2"   = "aCD3_aCD28_a4_1BB",
                                    "aCD3_aCD28_a4_1BB.3"   = "aCD3_aCD28_a4_1BB",
                                    "aCD3_aCD28_a4_1BB.4"   = "aCD3_aCD28_a4_1BB",
                                    "aCD3_aCD27_aCD28_a4_1BB.1"   = "aCD3_aCD27_aCD28_a4_1BB",
                                    "aCD3_aCD27_aCD28_a4_1BB.2"   = "aCD3_aCD27_aCD28_a4_1BB",
                                    "aCD3_aCD27_aCD28_a4_1BB.3"   = "aCD3_aCD27_aCD28_a4_1BB",
                                    "aCD3_aCD27_aCD28_a4_1BB.4"   = "aCD3_aCD27_aCD28_a4_1BB",
                                    "aCD3_aCD27_a4_1BB.1"   = "aCD3_aCD27_a4_1BB",
                                    "aCD3_aCD27_a4_1BB.2"   = "aCD3_aCD27_a4_1BB",
                                    "aCD3_aCD27_a4_1BB.3"   = "aCD3_aCD27_a4_1BB",
                                    "aCD3_aCD27_a4_1BB.4"   = "aCD3_aCD27_a4_1BB"
                                    # Voeg hier meer mapping toe zoals jij dat nodig hebt
  ))


#delete condtin t 

t_long <- t_long %>%
  select(-Condition_t)


# 3. Proteomics naar long + filter op genen
p_long <- p_z_score %>%
  filter(X %in% genes_of_interest) %>%
  pivot_longer(
    cols = -X,
    names_to = "Condition_p",
    values_to = "Proteomics_z"
  )

p_long <- p_long %>%
  mutate(condition_cleaned = recode(Condition_p,
                                    "X1.1.unstimulated" = "unstim",
                                    "X1.2.unstimulated" = "unstim",
                                    "X1.3.unstimulated" = "unstim",
                                    "X2.1.aCD3"   = "aCD3",
                                    "X2.2.aCD3"   = "aCD3",
                                    "X2.3.aCD3"   = "aCD3",
                                    "X3.1.aCD3.aCD27"   = "aCD3_aCD27",
                                    "X3.2.aCD3.aCD27"   = "aCD3_aCD27",
                                    "X3.3.aCD3.aCD27"   = "aCD3_aCD27",
                                    "X4.1.aCD3.a4_1BB"   = "aCD3_a4_1BB",
                                    "X4.2.aCD3.a4_1BB"   = "aCD3_a4_1BB",
                                    "X4.3.aCD3.a4_1BB"   = "aCD3_a4_1BB",
                                    "X5.2.aCD3.aCD28"   = "aCD3_aCD28",
                                    "X5.3.aCD3.aCD28"   = "aCD3_aCD28",
                                    "X5.1.aCD3.aCD28"   = "aCD3_aCD28",
                                    "X6.1.aCD3.aCD27.aCD28"   = "aCD3_aCD27_aCD28",
                                    "X6.2.aCD3.aCD27.aCD28"   = "aCD3_aCD27_aCD28",
                                    "X6.3.aCD3.aCD27.aCD28"   = "aCD3_aCD27_aCD28",
                                    "X7.1.aCD3.aCD28.a4_1BB"   = "aCD3_aCD28_a4_1BB",
                                    "X7.2.aCD3.aCD28.a4_1BB"   = "aCD3_aCD28_a4_1BB",
                                    "X7.3.aCD3.aCD28.a4_1BB"   = "aCD3_aCD28_a4_1BB"
                                    # Voeg hier meer mapping toe zoals jij dat nodig hebt
  ))


p_long <- p_long %>%
  select(-Condition_p)



#calculate the mean 

t_long_mean <- t_long %>% 
  group_by(X, condition_cleaned) %>% 
  summarise(mean_value_t = mean(Transcriptomics_z), .groups = "drop")



p_long_mean <- p_long %>% 
  group_by(X, condition_cleaned) %>% 
  summarise(mean_value_p = mean(Proteomics_z), .groups = "drop")



df_merged <- t_long_mean %>% 
  full_join(p_long_mean, by = c("X", "condition_cleaned"))

df_filtered <- df_merged %>%
  filter(condition_cleaned %in% c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB"))


df_filtered$condition_cleaned <- factor(df_filtered$condition_cleaned, levels = c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB"))



# Maak plots per conditie
conditions <- unique(df_filtered$condition_cleaned)

x_limits <- range(df_filtered$mean_value_t, na.rm = TRUE)
y_limits <- range(df_filtered$mean_value_p, na.rm = TRUE)

for (cond in conditions) {
  df_subset <- df_filtered %>% filter(condition_cleaned == cond)
  
  p <- ggplot(df_subset, aes(x = mean_value_t, y = mean_value_p)) +
    geom_point(size = 1, alpha = 0.7) +
    theme_classic() +
    xlim(x_limits) +
    ylim(y_limits) + 
    labs(
      title = paste("Scatterplot:", cond),
      x = "Mean transcriptomics z-score",
      y = "Mean proteomics z-score"
    )
  
  
  output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/density/"
  output_base <- paste0(output_dir, "new_DE_integration_", cond)
  
  # Sla op als PNG (hoge resolutie)
  ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)
  
  
  # Sla op als SVG
  ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)
  
  
  # Sla op als EPS (voor vectorformaat)
  ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)
  
  
}


 







### old code ### 

output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/density/"
output_base <- paste0(output_dir, "DE_integration")

# Sla op als PNG (hoge resolutie)
ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)


# Sla op als SVG
ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)


# Sla op als EPS (voor vectorformaat)
ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)



df_avg_cond <- df_filtered %>%
  group_by(X, condition_cleaned) %>%
  summarise(
    mean_Transcriptomics_z = mean(Transcriptomics_z, na.rm = TRUE),
    mean_Proteomics_z      = mean(Proteomics_z, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(df_avg_cond, aes(x = mean_Transcriptomics_z, y = mean_Proteomics_z)) +
  geom_point() +
  facet_wrap(~ condition_cleaned) +
  theme_minimal()



# 4. Join transcriptomics en proteomics op gen + uniforme conditienaam
df_merged <- t_long %>% 
  full_join(p_long, by = c("X", "condition_cleaned"))

df_filtered <- df_merged %>%
  filter(condition_cleaned %in% c("aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB"))



# 5. Maak een 2D-density-plot per conditie (facet)
ggplot(df_filtered, aes(x = Transcriptomics_z, y = Proteomics_z)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.4) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(~ condition_cleaned) +
  theme_minimal() +
  labs(
    title = "2D density plot: Transcriptomics vs Proteomics",
    x = "Transcriptomics z-score",
    y = "Proteomics z-score"
  )



# df_avg <- df_filtered %>%
#   group_by(X) %>%  # groepeer op gen
#   summarise(
#     mean_Transcriptomics_z = mean(Transcriptomics_z, na.rm = TRUE),
#     mean_Proteomics_z      = mean(Proteomics_z, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# 
# 
# library(ggplot2)
# 
# ggplot(df_avg, aes(x = mean_Transcriptomics_z, y = mean_Proteomics_z)) +
#   geom_point() +
#   theme_minimal() +
#   facet_wrap(~ condition_cleaned) +
#   labs(
#     x = "Mean Transcriptomics z-score",
#     y = "Mean Proteomics z-score"
#   )




df_avg_cond <- df_filtered %>%
  group_by(X, condition_cleaned) %>%
  summarise(
    mean_Transcriptomics_z = mean(Transcriptomics_z, na.rm = TRUE),
    mean_Proteomics_z      = mean(Proteomics_z, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(df_avg_cond, aes(x = mean_Transcriptomics_z, y = mean_Proteomics_z)) +
  geom_point() +
  facet_wrap(~ condition_cleaned) +
  theme_minimal()
