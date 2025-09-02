library(dplyr)
Pearson_transcript_gene <- read.csv("PearsonSpearman_CCLE_clean/Pearson_transcript_clean.csv")
Spearman_transcript_gene <- read.csv("PearsonSpearman_CCLE_clean/Spearman_transcript_clean.csv")

transcript_Pearson_20pct_higher <- Pearson_transcript_gene %>%
  filter(
    !is.na(Pearson_transcript) & !is.na(Pearson_gene),           # remove NAs
    Pearson_transcript >= 0.5,                                    # new condition
    abs(Pearson_gene) > 0,                                        # avoid divide-by-zero
    (Pearson_transcript - Pearson_gene) / abs(Pearson_gene) >= 0.2
  ) %>%
  pull(transcript_ID)

transcript_Spearman_20pct_higher <- Spearman_transcript_gene %>%
  filter(
    !is.na(Spearman_transcript) & !is.na(Spearman_gene),           # remove NAs
    Spearman_transcript >= 0.5,                                    # new condition
    abs(Spearman_gene) > 0,                                        # avoid divide-by-zero
    (Spearman_transcript - Spearman_gene) / abs(Spearman_gene) >= 0.2
  ) %>%
  pull(transcript_ID)

transcript_union <- unique(c(transcript_Spearman_20pct_higher, transcript_Pearson_20pct_higher))

# Check how many unique transcripts are included
length(transcript_union)

# View first few
head(transcript_union)

write.csv(transcript_union, "transcript_union_20pct_higher.csv", row.names = FALSE)

old_discordant <- read.csv("NEW_discordant_protein_coding.csv")
# Get intersecting transcripts
shared_transcripts <- intersect(old_discordant$transcript_ID, transcript_union)

# Count how many overlap
length(shared_transcripts)

# View first few overlaps
head(shared_transcripts)

transcript_union_annotated <- Pearson_transcript_gene %>%
  filter(transcript_ID %in% transcript_union)
shared_transcripts_annotated <- Pearson_transcript_gene %>%
  filter(transcript_ID %in% shared_transcripts)
write.csv(transcript_union_annotated, "transcript_union_annotated.csv", row.names = FALSE)
write.csv(shared_transcripts_annotated, "shared_transcripts_annotated.csv", row.names = FALSE)


