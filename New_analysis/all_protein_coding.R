library(readr)
library(dplyr)
library(ggplot2)

CCLE_trans_expr <- read_tsv("Data/CCLE_Melanoma_transcript_expression.txt")
CCLE_trans_expr <- CCLE_trans_expr %>%
  mutate(transcript_id = sub("\\..*$", "", transcript_id))
CCLE_trans_expr_clean <- CCLE_trans_expr[,-1]

# Keep only transcripts with >10 counts in at least 3 samples
filtered_transcripts <- CCLE_trans_expr_clean[
  rowSums(CCLE_trans_expr_clean[, 2:64] > 1) >= 2,
]

protein_coding <-read.csv("Data/NEW_all_protein_coding.csv")

# First, ensure transcript_id is character and strip version numbers (if needed)
filtered_transcripts$transcript_id <- sub("\\..*$", "", as.character(filtered_transcripts$transcript_id))

# If protein_coding is a data frame:
protein_coding$transcript_id <- sub("\\..*$", "", as.character(protein_coding$transcript_id))

# If protein_coding is a vector:
# protein_coding <- sub("\\..*$", "", as.character(protein_coding))

# Now filter
final_transcripts <- filtered_transcripts[
  filtered_transcripts$transcript_id %in% protein_coding$transcript_id,  # or just protein_coding if a vector
]

correlation <- read.csv("Data/NEW_correlation_protein_coding.csv")
final <- read.csv("Data/NEW_discordant_protein_coding.csv")

write.csv(final_transcripts,"Data/ALL_protein_coding_FINAL.csv")
