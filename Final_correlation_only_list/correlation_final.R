correlations <- read.csv("Data/correlation_list.csv")

# Assuming your data.frame is named 'correlations'

# Filter based on the criteria
filtered_transcripts <- correlations[
  (correlations$spearman_TCGA_value > 0.3 | correlations$pearson_TCGA_value > 0.3) &
    (correlations$spearman_CCLE_value > 0.5 | correlations$pearson_CCLE_value > 0.5),
]

write.csv(filtered_transcripts,"0.3_0.5_new_draft_list.csv")

# Extract the unique transcript IDs (assumes IDs are consistent across columns)
transcript_ids <- unique(c(
  filtered_transcripts$transcriptid_spearman_TCGA,
  filtered_transcripts$transcriptid_pearson_TCGA,
  filtered_transcripts$transcriptid_spearman_CCLE,
  filtered_transcripts$transcriptid_pearson_CCLE
))

# Clean NAs if any
transcript_ids <- transcript_ids[!is.na(transcript_ids)]

# Export to CSV
write.csv(data.frame(transcript_id = transcript_ids), "correlation_CCLE_0.5_TCGA_0.3_FINAL.csv", row.names = FALSE)


### stringent correlation
correlations <- read.csv("Data/correlation_list.csv")

# Assuming your data.frame is named 'correlations'

# Filter based on the criteria
filtered_transcripts_stringent <- correlations[
  (correlations$spearman_TCGA_value > 0.3 & correlations$pearson_TCGA_value > 0.3) &
    (correlations$spearman_CCLE_value > 0.5 & correlations$pearson_CCLE_value > 0.5),
]

# Extract the unique transcript IDs (assumes IDs are consistent across columns)
transcript_ids_stringent <- unique(c(
  filtered_transcripts_stringent$transcriptid_spearman_TCGA,
  filtered_transcripts_stringent$transcriptid_pearson_TCGA,
  filtered_transcripts_stringent$transcriptid_spearman_CCLE,
  filtered_transcripts_stringent$transcriptid_pearson_CCLE
))

# Clean NAs if any
transcript_ids_stringent <- transcript_ids_stringent[!is.na(transcript_ids_stringent)]

# Export to CSV
write.csv(data.frame(transcript_ids_stringent), "correlation_CCLE_0.5_TCGA_0.3_stringent_FINAL.csv", row.names = FALSE)

