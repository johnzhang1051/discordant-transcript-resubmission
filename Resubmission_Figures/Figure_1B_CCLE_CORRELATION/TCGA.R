library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

expr_df <- read.csv("Data/Transcript_expression_melanoma_log2.csv")
final_discordant <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")





### for paper plot
ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.1) +
  scale_color_manual(
    values = c("transcript_unique" = "blue", "transcript_all_other" = "red")
  ) +
  labs(
    title = "ABR vs MITF Expression",
    x = "MITF expression (log2 + 1)",
    y = "ABR expression (log2 + 1)",
    color = "Transcript group"
  ) +
  theme_classic(base_size = 13, base_family = "Times") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
