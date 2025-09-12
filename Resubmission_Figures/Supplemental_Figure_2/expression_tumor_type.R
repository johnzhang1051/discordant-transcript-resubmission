library(ggplot2)
library(dplyr)
library(forcats)

#comprehensive_transcript_dataset <- read.csv("data/comprehensive_transcript_dataset.csv")
#comprehensive_transcript_dataset <- comprehensive_transcript_dataset %>% filter(is_protein_coding_filtered)

### ABR Gene expression CCLE - Top 20 tumor types
CCLE_gene_for_figure <- readRDS("data/CCLE_gene_disease_Figure1B_new.rds")
colnames(CCLE_gene_for_figure)[3:ncol(CCLE_gene_for_figure)] <- gsub(" \\(\\d+\\)$", "", colnames(CCLE_gene_for_figure)[3:ncol(CCLE_gene_for_figure)])

# Prepare the data and get top 20 tumor types by median ABR
abr_df <- CCLE_gene_for_figure %>%
  dplyr::select(DepmapModelType = 3, ABR) %>%
  group_by(DepmapModelType) %>%
  summarise(median_ABR = median(ABR, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_ABR)) %>%
  #slice_head(n = 20) %>%
  pull(DepmapModelType)

# Filter to top 20 and prepare for plotting
abr_df_plot <- CCLE_gene_for_figure %>%
  dplyr::select(DepmapModelType = 3, ABR) %>%
  filter(DepmapModelType %in% abr_df) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df_plot <- abr_df_plot %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ABR, .fun = median))

# Plot
ggplot(abr_df_plot, aes(x = DepmapModelType, y = ABR, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "blue", "Other" = "gray70")) +
  labs(
    title = "ABR Expression by Tumor Type",
    x = "DepMap Model Type",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

### ABR TRANSCRIPT EXPRESSION CCLE - ENST00000544583
### import df and plot ABR
filtered_df_for_ABR_CCLE_transcript <- readRDS("data/Figure_1C_CCLE_transcript_new.rds")

# Get top 20 tumor types by median transcript expression
abr_transcript_top20 <- filtered_df_for_ABR_CCLE_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000544583) %>%
  group_by(DepmapModelType) %>%
  summarise(median_transcript = median(ENST00000544583, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_transcript)) %>%
  slice_head(n = 20) %>%
  pull(DepmapModelType)

# Prepare the data
abr_transcript_plot <- filtered_df_for_ABR_CCLE_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000544583) %>%
  filter(DepmapModelType %in% abr_transcript_top20) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_transcript_plot <- abr_transcript_plot %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000544583, .fun = median))

# Plot
ggplot(abr_transcript_plot, aes(x = DepmapModelType, y = ENST00000544583, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Transcript Expression Across Top 20 Tumor Types (by median expression)",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Transcript Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

################## NEW PLOTS


### PEX10 TRANSCRIPT EXPRESSION CCLE - ENST00000508384
### import df and plot PEX10
filtered_df_for_transcript <- readRDS("data/Figure_1C_CCLE_transcript_new.rds")

# Get top 20 tumor types by median PEX10 transcript expression
pex10_transcript_top20 <- filtered_df_for_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000508384) %>%
  group_by(DepmapModelType) %>%
  summarise(median_transcript = median(ENST00000508384, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_transcript)) %>%
  #slice_head(n = 20) %>%
  pull(DepmapModelType)

# Prepare the data for PEX10
pex10_transcript_plot <- filtered_df_for_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000508384) %>%
  filter(DepmapModelType %in% pex10_transcript_top20) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median PEX10 (highest on right)
pex10_transcript_plot <- pex10_transcript_plot %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000508384, .fun = median))

# Plot PEX10
ggplot(pex10_transcript_plot, aes(x = DepmapModelType, y = ENST00000508384, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "PEX10 Transcript Expression Across Top 20 Tumor Types (by median expression)",
    x = "DepMap Model Type (ordered by median PEX10 expression)",
    y = "PEX10 Transcript Expression (ENST00000508384)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


### METTL9 TRANSCRIPT EXPRESSION CCLE - ENST00000562961
# Get top 20 tumor types by median METTL9 transcript expression
mettl9_transcript_top20 <- filtered_df_for_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000562961) %>%
  group_by(DepmapModelType) %>%
  summarise(median_transcript = median(ENST00000562961, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_transcript)) %>%
  #slice_head(n = 20) %>%
  pull(DepmapModelType)

# Prepare the data for METTL9
mettl9_transcript_plot <- filtered_df_for_transcript %>%
  dplyr::select(DepmapModelType = 3, ENST00000562961) %>%
  filter(DepmapModelType %in% mettl9_transcript_top20) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median METTL9 (highest on right)
mettl9_transcript_plot <- mettl9_transcript_plot %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000562961, .fun = median))

# Plot METTL9
ggplot(mettl9_transcript_plot, aes(x = DepmapModelType, y = ENST00000562961, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "METTL9 Transcript Expression Across Top 20 Tumor Types (by median expression)",
    x = "DepMap Model Type (ordered by median METTL9 expression)",
    y = "METTL9 Transcript Expression (ENST00000562961)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# Optional: Check if SKCM is in top 20 for each transcript
print(paste("SKCM included in PEX10 plot:", "SKCM" %in% pex10_transcript_top20))
print(paste("SKCM included in METTL9 plot:", "SKCM" %in% mettl9_transcript_to
