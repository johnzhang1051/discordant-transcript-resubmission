library(ggplot2)
library(dplyr)

# Load main data
Gene_transcript_pearson <- read.csv("Data/merged_pearsons_annotated.csv")
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
ggplot(Gene_transcript_pearson, aes(x = pearson_gene, y = pearson_transcript)) +
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
