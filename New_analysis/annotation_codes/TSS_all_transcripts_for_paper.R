# Load required package
library(tidyverse)
final_transcript <-read.csv("Data/final_transcript_list.txt")
correlation_only <- read.csv("Data/correlation_only_list.txt")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")
# Assume your input dataframe is named `tss_data`
# It should already contain: gene_id, transcript_id, TSS, Strand, Chromosome


# Ensure data types are consistent
tss_data <- tss_data %>%
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  )

# Step 1: Sort TSS by chromosome and TSS position
tss_with_spacing <- tss_data %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),       # Distance from previous TSS
    next_dist = abs(lead(TSS) - TSS),      # Distance to next TSS
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)  # Closest neighbor
  ) %>%
  ungroup()

# Step 2: Label whether TSS is isolated (≥ 500bp from any other)
tss_with_spacing <- tss_with_spacing %>%
  mutate(is_isolated = min_dist >= 500)

# Step 3: Select only transcript and isolation status
result <- tss_with_spacing %>%
  select(transcript_id, is_isolated)

# Optional: View or export the result
print(result)
# write_tsv(result, "Results/transcript_isolation_flags.tsv")


# Count how many TRUE and FALSE values
summary_counts <- result %>%
  count(is_isolated)

# Calculate percentages
summary_percent <- summary_counts %>%
  mutate(percent = 100 * n / sum(n))

# View the summary
print(summary_percent)




### protein coding only
# Load required package
library(tidyverse)

# Load data
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")

# Merge with transcript_type and keep only protein_coding
tss_data_protein_coding <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Ensure data types are consistent
tss_data_protein_coding <- tss_data_protein_coding %>%
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  )

# Step 1: Sort TSS by chromosome and TSS position
tss_with_spacing <- tss_data_protein_coding %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),       # Distance from previous TSS
    next_dist = abs(lead(TSS) - TSS),      # Distance to next TSS
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)  # Closest neighbor
  ) %>%
  ungroup()

# Step 2: Label whether TSS is isolated (≥ 50bp from any other)
tss_with_spacing <- tss_with_spacing %>%
  mutate(is_isolated = min_dist >= 500)

# Step 3: Select only transcript and isolation status
result <- tss_with_spacing %>%
  select(transcript_id, is_isolated)

# Optional: View or export the result
print(result)
# write_tsv(result, "Results/transcript_isolation_flags.tsv")

# Count how many TRUE and FALSE values
summary_counts <- result %>%
  count(is_isolated)

# Calculate percentages
summary_percent <- summary_counts %>%
  mutate(percent = 100 * n / sum(n))

# View the summary
print(summary_percent)

