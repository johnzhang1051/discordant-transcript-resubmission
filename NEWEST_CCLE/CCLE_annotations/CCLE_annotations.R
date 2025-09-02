library(data.table)
library(dplyr)

# Load the file
CCLE_melanoma_transcript <- fread("Transcript_expression_melanoma_log2.csv")
transcript_type <-read.csv("transcripttype.csv")
transcript_id <- (transcript_type %>%
  filter(transcript_type == "protein_coding") %>%
  pull(transcript_id))
transcript_id <- as.data.frame(transcript_id)

# Set first column name to 'sample'
setnames(CCLE_melanoma_transcript, 1, "sample")

# Compute column-wise mean (i.e., average log2 expression per transcript)
# Exclude the first column (sample names)
transcript_means <- colMeans(CCLE_melanoma_transcript[, -1, with = FALSE], na.rm = TRUE)

# Create final data.table with transcript_id and average log2 expression
CCLE_melanoma_transcript_averageexpression <- data.table(
  transcript_id = names(transcript_means),
  avg_log2_expression = as.numeric(transcript_means)
)


CCLE_melanoma_average_protein_coding <- merge(CCLE_melanoma_transcript_averageexpression, transcript_id, by="transcript_id")

CCLE_Pearson <- read.csv("Pearson_transcript_clean.csv")
CCLE_Spearman <- read.csv("Spearman_transcript_clean.csv")