
# Load required package
library(tidyverse)

# Step 1: Load transcript list (one transcript ID per line)
transcript_list <- read_tsv("Data/CCLE_final_merged.csv", col_names = "transcript_id")

# Step 2: Load and clean TSS data
tss_data <- read_tsv("Data/biomart_TSS_chrom_simple.txt") %>%
  rename_with(~ str_replace_all(., "\\s+", "_")) %>%  # Clean up any leftover spaces
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  ) %>%
  select(gene_id, transcript_id, TSS, Strand, Chromosome)

# Step 3: Compute spacing between TSSs on each chromosome
tss_with_spacing <- tss_data %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),
    next_dist = abs(lead(TSS) - TSS),
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 4: Tag isolated TSS (TRUE if ≥500bp away from neighbors)
tss_with_spacing <- tss_with_spacing %>%
  mutate(is_isolated = min_dist >= 500)

# Step 5: For each transcript in the input list, check if it’s isolated
result <- transcript_list %>%
  left_join(tss_with_spacing %>% select(transcript_id, is_isolated), by = "transcript_id") %>%
  mutate(is_isolated = replace_na(is_isolated, FALSE))  # transcripts not found → FALSE

# Step 6: Write results to file
write_csv(result, "transcript_distal.csv")

# Optional: view results
print(result)
