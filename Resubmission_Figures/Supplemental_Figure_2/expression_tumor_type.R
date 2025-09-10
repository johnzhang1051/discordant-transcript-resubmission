library(ggplot2)
library(dplyr)
comprehensive_transcript_dataset <- read.csv("data/comprehensive_transcript_dataset.csv")

### ABR Gene expression CCLE

library(ggplot2)
library(dplyr)
library(forcats)  # for fct_reorder
## 
CCLE_gene_for_figure <- readRDS("data/CCLE_gene_disease_Figure1B_new.rds")
colnames(CCLE_gene_for_figure)[3:ncol(CCLE_gene_for_figure)] <- gsub(" \\(\\d+\\)$", "", colnames(CCLE_gene_for_figure)[3:ncol(CCLE_gene_for_figure)])
# Prepare the data
abr_df <- CCLE_gene_for_figure %>%
  dplyr::select(DepmapModelType = 3, ABR) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ABR, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ABR, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "blue", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



### ABR TRANSCRIPT EXPRESSION CCLE
### import df and plot ABR
filtered_df_for_ABR_CCLE_transcript <- readRDS("data/Figure_1C_CCLE_transcript_new.rds")


# Prepare the data
abr_df <- filtered_df_for_ABR_CCLE_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000544583) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000544583, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ENST00000544583, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
