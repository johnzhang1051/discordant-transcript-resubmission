# Load required package
library(tidyverse)

# Step 0: Load transcript type data and TSS data
transcript_type <- read_csv("Data/transcripttype.csv", col_names = c("transcript_id", "transcript_type"))
tss_data <- read_csv("Data/biomart_TSS_chrom_simple.csv")

# Step 1: Join TSS data with transcript type and filter for protein_coding only
tss_data <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Step 2: Ensure data types are consistent
tss_data <- tss_data %>%
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  )

# Step 3: Sort TSS by chromosome and TSS position, compute spacing
tss_with_spacing <- tss_data %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),       # Distance from previous TSS
    next_dist = abs(lead(TSS) - TSS),      # Distance to next TSS
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)  # Closest neighbor
  ) %>%
  ungroup()

# Step 4: Label whether TSS is isolated (≥ 500bp from any other)
tss_with_spacing <- tss_with_spacing %>%
  mutate(is_isolated = min_dist >= 500)

# Step 5: Select only transcript and isolation status
result <- tss_with_spacing %>%
  select(transcript_id, is_isolated)

# Step 6: Print results
print(result)

# Step 7: Count how many TRUE and FALSE values
summary_counts <- result %>%
  count(is_isolated)

# Step 8: Calculate percentages
summary_percent <- summary_counts %>%
  mutate(percent = 100 * n / sum(n))

# Step 9: View the summary
print(summary_percent)

# Export the result to a TSV file
write_tsv(result, "distal_TSS_protein_coding.tsv")


