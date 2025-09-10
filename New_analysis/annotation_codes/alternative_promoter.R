library(tidyverse)

# Load data
final_transcript <-read.csv("Data/Discordant_RESUBMISSION.csv")
correlation_transcript <- read.csv("Data/correlated_RESUBMISSION.csv")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")

promoter_by_transcript <- read.csv("Data_promoter/promoter_by_transcript.csv")
# Remove the period and following numbers from transcript_id
promoter_by_transcript$transcript_id <- sub("\\..*$", "", promoter_by_transcript[[1]])
promoter_by_transcript <- promoter_by_transcript %>%
  group_by(geneId, promoterId) %>%
  mutate(promoter_count_within_gene = n()) %>%  # how many times this promoter appears in same gene
  ungroup() %>%
  group_by(geneId, transcript_id) %>%
  mutate(is_unique_promoter = promoter_count_within_gene == 1) %>%
  ungroup()
promoter_by_transcript <- promoter_by_transcript[,-(2:4)]


# Merge final transcripts
final_annotated <- final_transcript %>%
  inner_join(promoter_by_transcript, by = "transcript_id")

# Merge correlation-only transcripts
correlation_annotated <- correlation_transcript %>%
  inner_join(promoter_by_transcript, by = "transcript_id")


library(dplyr)
library(ggplot2)

# Step 1: Add group labels
all_df <- promoter_by_transcript %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "All Transcripts")

correlation_df <- correlation_annotated %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "Correlation Only")

final_df <- final_annotated %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "Final Transcript")

# Step 2: Combine all into one long dataframe
combined_df <- bind_rows(all_df, correlation_df, final_df)

# Step 3: Summarize for plotting
summary_df <- combined_df %>%
  group_by(group) %>%
  summarise(
    proportion_unique = mean(is_unique_promoter, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Plot
ggplot(summary_df, aes(x = group, y = proportion_unique, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(proportion_unique, 2)), vjust = -0.5, size = 4) +
  ylim(0, 1.1) +
  labs(
    title = "Proportion of Transcripts with Unique Promoters",
    x = "Transcript Group",
    y = "Proportion Unique"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



### protein coding transcripts
library(tidyverse)

# Load data
final_transcript <-read.csv("Data/Discordant_RESUBMISSION.csv")
correlation_transcript <- read.csv("Data/correlated_RESUBMISSION.csv")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")
correlation_transcript <- correlation_transcript %>%
  filter(!(transcript_id %in% final_transcript$transcript_id))

promoter_by_transcript <- read.csv("Data_promoter/promoter_by_transcript.csv")

# Clean transcript IDs
promoter_by_transcript$transcript_id <- sub("\\..*$", "", promoter_by_transcript$transcript_id)
transcript_type$transcript_id <- sub("\\..*$", "", transcript_type$transcript_id)

# Keep only protein_coding transcripts
protein_coding_ids <- transcript_type %>%
  filter(transcript_type == "protein_coding") %>%
  pull(transcript_id)

promoter_by_transcript <- promoter_by_transcript %>%
  filter(transcript_id %in% protein_coding_ids) %>%
  group_by(geneId, promoterId) %>%
  mutate(promoter_count_within_gene = n()) %>%
  ungroup() %>%
  group_by(geneId, transcript_id) %>%
  mutate(is_unique_promoter = promoter_count_within_gene == 1) %>%
  ungroup() %>%
  select(transcript_id, is_unique_promoter)

# Merge final and correlation transcript sets (filtering to protein_coding only)
final_annotated <- final_transcript %>%
  filter(transcript_id %in% protein_coding_ids) %>%
  inner_join(promoter_by_transcript, by = "transcript_id")

correlation_annotated <- correlation_transcript %>%
  filter(transcript_id %in% protein_coding_ids) %>%
  inner_join(promoter_by_transcript, by = "transcript_id")

# Prepare all three groups
all_df <- promoter_by_transcript %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "All Transcripts")

correlation_df <- correlation_annotated %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "Correlation Only")

final_df <- final_annotated %>%
  select(transcript_id, is_unique_promoter) %>%
  mutate(group = "Final Transcript")

# Combine all into one long dataframe
combined_df <- bind_rows(all_df, correlation_df, final_df)

# Summarize for plotting
summary_df <- combined_df %>%
  group_by(group) %>%
  summarise(
    proportion_unique = mean(is_unique_promoter, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(summary_df, aes(x = group, y = proportion_unique, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(proportion_unique, 2)), vjust = -0.5, size = 4) +
  ylim(0, 1.1) +
  labs(
    title = "Proportion of Protein-Coding Transcripts with Unique Promoters",
    x = "Transcript Group",
    y = "Proportion Unique"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Plot with custom colors
# Final plot with bars closer together and no labels
ggplot(summary_df, aes(x = group, y = proportion_unique, fill = group)) +
  geom_bar(stat = "identity", width = 0.8) +  # Wider bars
  ylim(0, 1.1) +
  labs(
    title = "Proportion of Protein-Coding Transcripts with Unique Promoters",
    x = "Transcript Group",
    y = "Proportion Unique"
  ) +
  scale_fill_manual(values = c(
    "All Transcripts" = "#FADADD",    # light pink
    "Correlation Only" = "#F08080",   # medium pink
    "Final Transcript" = "#FF0000"    # red
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

make_fisher_table <- function(group1, group2, data) {
  data_filtered <- data %>%
    filter(group %in% c(group1, group2))
  
  counts <- table(data_filtered$group, data_filtered$is_unique_promoter)
  
  # Ensure correct row and column ordering
  result <- matrix(
    c(counts[group1, "TRUE"],  counts[group1, "FALSE"],
      counts[group2, "TRUE"],  counts[group2, "FALSE"]),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Group = c(group1, group2),
      UniquePromoter = c("TRUE", "FALSE")
    )
  )
  
  return(result)
}
