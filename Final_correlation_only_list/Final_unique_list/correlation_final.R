## THIS CODE GENERATES 

library(dplyr)

CCLE_pearson <- read.csv("Data/CCLE_pearson.csv")
CCLE_spearman <- read.csv("Data/CCLE_spearman.csv")
TCGA_pearson <- read.csv("Data/TCGA_pearson.csv")
TCGA_spearman <- read.csv("Data/TCGA_spearman.csv")


# Filter CCLE_pearson for desired conditions
filtered_CCLE_pearson <- CCLE_pearson %>%
  filter(
    coefficient_transcript >= 0.5,
    coefficient_transcript >= 1.2 * coefficient_gene
  )

# View or save
head(filtered_CCLE)

filtered_CCLE_spearman <- CCLE_spearman %>%
  filter(
    coefficient_transcript >= 0.5,
    coefficient_transcript >= 1.2 * coefficient_gene
  )

merged_CCLE <- full_join(
  filtered_CCLE_spearman,
  filtered_CCLE_pearson,
  by = "transcript_id"
)

# Filter each TCGA list
filtered_TCGA_pearson <- TCGA_pearson %>%
  filter(coefficient_transcript >= 0.3)

filtered_TCGA_spearman <- TCGA_spearman %>%
  filter(coefficient_transcript >= 0.3)

# Merge to keep transcripts found on either list
merged_TCGA <- full_join(
  filtered_TCGA_pearson,
  filtered_TCGA_spearman,
  by = "transcript_id",
  suffix = c("_pearson", "_spearman")
)


#### Get the shared transcript_ids
final_transcripts <- as.data.frame(intersect(merged_CCLE$transcript_id, merged_TCGA$transcript_id))
colnames(final_transcripts) <- "transcript_id"


write.csv(final_transcripts, "Final_CCLE_TCGA_overlap.csv", row.names = FALSE)

transcript_type <- read.csv("Data/transcripttype.csv")

### protein coding of discorantly correlated
# First, merge transcript_type info into final_transcript
final_transcript_merged <- merge(
  final_transcripts,
  transcript_type,
  by = "transcript_id"
)

# Then, filter to keep only protein_coding transcripts
final_transcript_protein_coding <- final_transcript_merged %>%
  filter(transcript_type == "protein_coding")


##Protein coding for all transcripts
# Subset transcript_type to keep only protein_coding entries
all_transcript_protein_coding <- transcript_type %>%
  filter(transcript_type == "protein_coding")

## protein coding of correlated list
correlation_only <- read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")
correlation_transcript_merged <- merge(
  correlation_only,
  transcript_type,
  by = "transcript_id"
)

# Then, filter to keep only protein_coding transcripts
correlation_transcript_protein_coding <- correlation_transcript_merged %>%
  filter(transcript_type == "protein_coding")

write.csv(all_transcript_protein_coding,"NEW_all_protein_coding.csv")
write.csv(final_transcript_protein_coding,"NEW_discordant_protein_coding.csv")
write.csv(correlation_transcript_protein_coding,"NEW_correlation_protein_coding.csv")

