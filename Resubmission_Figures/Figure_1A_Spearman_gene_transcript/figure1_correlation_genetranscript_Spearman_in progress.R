library(ggplot2)
library(dplyr)

# Load main data
Gene_transcript_spearman <- read.csv("Data/Spearman_transcript_clean.csv")


# Load list of transcript IDs to highlight in pink
pink_transcripts_df <- read.csv("Data/Spearman_0.5.csv")
Tsoi_CCLE_correlated <- read.csv("Data/MITF_correlated_CCLE__Tsoi.csv")
pink_transcripts_df <-merge(pink_transcripts_df, Tsoi_CCLE_correlated, by="transcript_ID")

# Load list of transcript IDs to highlight in red (replacing the old hardcoded list)
red_transcripts_df <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")
Spearman_all <- read.csv("Data/Spearman_discordant.csv")
red_transcripts_df <- merge(red_transcripts_df, Spearman_all, by = "transcript_ID")

# Annotate highlight status
# Annotate highlight status (corrected)
Gene_transcript_spearman <- Gene_transcript_spearman %>%
  mutate(highlight = case_when(
    transcript_ID %in% red_transcripts_df$transcript_ID ~ "Red",
    transcript_ID %in% pink_transcripts_df$transcript_ID ~ "Pink",
    TRUE ~ "Other"
  ))


# Plot
ggplot(Gene_transcript_spearman, aes(x = Spearman_gene, y = Spearman_transcript)) +
  # Background points
  geom_point(
    data = Gene_transcript_spearman %>% filter(highlight == "Other"),
    color = "gray70", size = 0.1, alpha = 0.3
  ) +
  # Pink points
  geom_point(
    data = Gene_transcript_spearman %>% filter(highlight == "Pink"),
    color = "salmon", size = 0.1
  ) +
  # Red points
  geom_point(
    data = Gene_transcript_spearman %>% filter(highlight == "Red"),
    color = "red", size = 1.5
  ) +
  labs(
    title = "Gene vs Transcript Spearman Correlation",
    x = "Gene Spearman Correlation",
    y = "Transcript Spearman Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),        # Axis tick labels
    axis.title = element_text(size = 16),       # Axis titles
    plot.title = element_text(size = 18, hjust = 0.5)  # Plot title
)