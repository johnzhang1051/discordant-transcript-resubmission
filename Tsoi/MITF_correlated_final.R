### transcripts with 0.5 correlation in TSOI
Tsoi_0.5 <- read.csv("Data_final_MITF_Correlation/correlated_transcripts_r0.5.csv")
colnames(Tsoi_0.5)[1:3] <- c("transcript_ID", "Pearson_Tsoi", "Spearman_Tsoi")


### transcripts with 0.5 and CoCor significance in CCLE
CCLE_Spearman <- read.csv("Data_final_MITF_Correlation/CCLE_Spearman_transcript_clean.csv")
colnames(CCLE_Spearman)[2] <- ("Spearman_CCLE")
CCLE_Pearson <- read.csv("Data_final_MITF_Correlation/CCLE_Pearson_transcript_clean.csv")
colnames(CCLE_Pearson)[2] <- ("Pearson_CCLE")
### final MITF correlated list

# Step 1: Merge Tsoi_0.5 with CCLE_Pearson
tsoi_merged_pearson <- merge(Tsoi_0.5, CCLE_Pearson, by.x = "transcript_ID", by.y = "transcript_ID", all.x = TRUE)

# Step 2: Merge result with CCLE_Spearman
tsoi_full_merged <- merge(tsoi_merged_pearson, CCLE_Spearman, by.x = "transcript_ID", by.y = "transcript_ID", all.x = TRUE)

# Filter for CCLE_Pearson or CCLE_Spearman ≥ 0.5
final_filtered_0.5 <- tsoi_full_merged[
  tsoi_full_merged$Pearson_CCLE >= 0.5 | tsoi_full_merged$Spearman_CCLE >= 0.5,
]

# View result
head(final_filtered_0.5)

transcript_type <- read.csv("Data_final_MITF_Correlation/transcripttype.csv")

library(dplyr)

# Merge transcript_type into final_filtered_0.5
merged_df <- final_filtered_0.5 %>%
  left_join(transcript_type, by = "transcript_ID")

# Filter to keep only protein_coding transcripts
final_filtered_protein_coding <- merged_df %>%
  filter(transcript_type == "protein_coding")

# View result
head(final_filtered_protein_coding)

write.csv(final_filtered_protein_coding, "final_MITF_correlated.csv")

