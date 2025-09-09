library(ggplot2)
library(dplyr)

# Load main data
Gene_transcript_pearson <- read.csv("Data/Pearson_transcript_clean.csv")

# Load list of transcript IDs to highlight in pink
### needto have pearson specific and spearman specific form
pink_transcripts_df <- read.csv("Data/Pearson_0.5.csv")
Tsoi_CCLE_correlated <- read.csv("Data/MITF_correlated_CCLE__Tsoi.csv")
pink_transcripts_df <-merge(pink_transcripts_df, Tsoi_CCLE_correlated, by="transcript_ID")
write.csv(pink_transcripts_df, "MITF_correlated_exclude_unique.csv") ### for final MITF correlated list

# Load list of transcript IDs to highlight in red (replacing the old hardcoded list)
red_transcripts_df <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")
Pearson_all <- read.csv("Data/Pearson_discordant.csv")
red_transcripts_df <- merge(red_transcripts_df, Pearson_all, by = "transcript_ID")

# Annotate highlight status
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    transcript_ID %in% red_transcripts_df$transcript_ID ~ "Red",
    transcript_ID %in% pink_transcripts_df$transcript_ID ~ "Pink",
    TRUE ~ "Other"
  ))

# Plot
ggplot(Gene_transcript_pearson, aes(x = Pearson_gene, y = Pearson_transcript)) +
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),        # Axis tick labels
    axis.title = element_text(size = 16),       # Axis titles
    plot.title = element_text(size = 18, hjust = 0.5)  # Plot title
  )
