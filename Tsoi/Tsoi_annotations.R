Tsoi_melanoma_transcript <- read.csv("kallisto_transcript_TPM.csv")
# Step 1: Make sure column 1 is named 'transcript_id'
colnames(Tsoi_melanoma_transcript)[1] <- "transcript_id"

# Step 2: Calculate the average TPM per transcript (excluding the transcript_id column)
Tsoi_melanoma_transcript$avg_TPM <- rowMeans(Tsoi_melanoma_transcript[ , -1], na.rm = TRUE)

# Step 3: Log2 transform with +2 pseudocount
Tsoi_melanoma_avg_exp <- Tsoi_melanoma_transcript %>%
  select(transcript_id, avg_TPM) %>%
  mutate(log2_avg_TPM = log2(avg_TPM + 2))

Tsoi_correlation <- read.csv("correlation_to_ENST00000394351.csv")
