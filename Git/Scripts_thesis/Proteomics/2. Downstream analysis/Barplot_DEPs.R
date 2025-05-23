#make bar_plot for the DEPs. 

library(readxl)
library(dplyr)

#read excell
df_diverging_prot <- read_xlsx("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/proteomics_data.xlsx")

# make counts numeric
df_diverging_prot$Count <- as.numeric(df_diverging_prot$Count)

# desired order


desired_order <- c("aCD3_aCD27 - aCD3", 
                   "aCD3_a4_1BB - aCD3", 
                   "aCD3_aCD28 - aCD3",
                   "aCD3_aCD27_aCD28 - aCD3", 
                   "aCD3_aCD28_a4_1BB - aCD3")

df_diverging_prot$Comparison <- factor(df_diverging_prot$Comparison, levels = desired_order)


library(ggplot2)

p <- ggplot(df_diverging_prot, aes(x = Count, y = Comparison, fill = Regulation)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(title = "DE proteins pvalue < 0.05, abd logfc 0.5 compared to αCD3",
       x = "Number of proteins",
       y = "Comparison") +
  theme_classic()

# plot
print(p)



#save files
output_dir <- "I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/figures/deg/DE_prot"
output_base <- paste0(output_dir, "DE_proteomics")

ggsave(paste0(output_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)
ggsave(paste0(output_base, ".svg"), plot = p, width = 8, height = 6)
ggsave(paste0(output_base, ".eps"), plot = p, width = 8, height = 6, device = cairo_ps)
