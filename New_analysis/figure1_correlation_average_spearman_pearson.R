pearson <- read.csv("Data/merged_pearsons_annotated.csv")
spearmans <-read.csv("Data/merged_spearman_annotated_stringent.csv")
# Merge data frames by transcript_id
merged_df <- merge(
  pearson[, c("transcript_id", "gene_id", "Gene_name", "pearson_transcript", "pearson_gene")],
  spearmans[, c("transcript_id", "spearman_transcript", "spearman_gene")],
  by = "transcript_id"
)

# Calculate average correlations
merged_df$avg_transcript_correlation <- rowMeans(merged_df[, c("pearson_transcript", "spearman_transcript")], na.rm = TRUE)
merged_df$avg_gene_correlation <- rowMeans(merged_df[, c("pearson_gene", "spearman_gene")], na.rm = TRUE)

# Keep only desired columns
final_df <- merged_df[, c("transcript_id", "gene_id", "Gene_name", "avg_transcript_correlation", "avg_gene_correlation")]

# View result
head(final_df)

library(ggplot2)
library(dplyr)

# Load main data
Gene_transcript_pearson <- final_df
Gene_transcript_pearson$transcript_id <- sub("\\..*$", "", Gene_transcript_pearson$transcript_id)
Gene_transcript_pearson$gene_id <- sub("\\..*$", "", Gene_transcript_pearson$gene_id)

# Load list of transcript IDs to highlight in pink
pink_transcripts_df <- read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")
pink_transcripts <- unique(pink_transcripts_df$transcript_id)
pink_transcripts <- sub("\\..*$", "", pink_transcripts)

# Load list of transcript IDs to highlight in red (replacing the old hardcoded list)
red_transcripts_df <- read.csv("Data/Final_CCLE_TCGA_overlap.csv")
red_transcripts <- unique(red_transcripts_df$transcript_id)
red_transcripts <- sub("\\..*$", "", red_transcripts)

# Annotate highlight status
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    transcript_id %in% red_transcripts ~ "Red",
    transcript_id %in% pink_transcripts ~ "Pink",
    TRUE ~ "Other"
  ))

# Plot
ggplot(Gene_transcript_pearson, aes(x = avg_gene_correlation, y = avg_transcript_correlation)) +
  # Background points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Other"),
    color = "gray70", size = 0.1, alpha = 0.3
  ) +
  # Pink points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Pink"),
    color = "salmon", size = 0.1
  ) +
  # Red points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Red"),
    color = "red", size = 1.5
  ) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()
