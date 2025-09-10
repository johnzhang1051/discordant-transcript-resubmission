# Load required libraries
library(tidyverse)

# Step 1: Load transcripts RDS
transcripts <- readRDS("Data/hg38_annotated_transcripts_for_analysis_FINAL.rds")

# Step 2: Load list of transcript IDs to keep
transcript_list <- read_tsv("Data/list.txt", col_names = "transcript_id")

# Step 3: Subset transcripts to only those in the list
transcripts <- transcripts %>%
  semi_join(transcript_list, by = "transcript_id")

# Step 4: Ensure required columns are numeric
transcripts <- transcripts %>%
  mutate(
    MITF_peak_n = as.numeric(MITF_peak_n),
    EBOX_n = as.numeric(EBOX_n)
  )

# Step 5: Subset protein_coding transcripts for one of the plots
transcripts_protein_coding <- transcripts %>%
  filter(transcript_type == "protein_coding")

# Step 6: Create output directory if needed
dir.create("Results", showWarnings = FALSE)

# Step 7: Write PDF with all plots
pdf("hg38_transcript_feature_plots_subset.pdf", width = 8, height = 6)

# Plot 1: MITF_peak_n
p1 <- ggplot(transcripts, aes(x = MITF_peak_n)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "Distribution of MITF Peaks per Transcript (Subset)",
       x = "MITF_peak_n", y = "Count") +
  theme_minimal()
print(p1)

# Plot 2: EBOX_n
p2 <- ggplot(transcripts, aes(x = EBOX_n)) +
  geom_histogram(binwidth = 1, fill = "darkgreen", color = "black") +
  scale_x_continuous(breaks = 0:9, limits = c(0, 9)) +
  labs(title = "Distribution of EBOX Count per Transcript (Subset)",
       x = "EBOX_n", y = "Count") +
  theme_minimal()
print(p2)


# Plot 4: Remote_TSS for protein_coding only
p3 <- ggplot(transcripts_protein_coding, aes(x = distal_protein_coding_TSS)) +
  geom_bar(fill = "orange", color = "black") +
  labs(title = "Remote TSS (Protein-Coding Only, Subset)",
       x = "Remote_TSS", y = "Count") +
  theme_minimal()
print(p3)

# Step 8: Close PDF
dev.off()



# Load required libraries
library(tidyverse)

# Step 1: Load data
transcripts <- readRDS("Data/hg38_annotated_transcripts_for_analysis_FINAL.rds")
transcript_list <- read_tsv("Data/list.txt", col_names = "transcript_id")

# Step 2: Subset to transcripts in list
transcripts <- transcripts %>%
  semi_join(transcript_list, by = "transcript_id")

# Step 3: Clean and convert expression ratio columns
expression_cols <- c("GSE61967.H3A.y", "PRJEB30337.y", "GSE61967.501mel.y", "Henja.y")

# Step 4: Replace "0" and "#DIV/0!" with NA, then convert to numeric
transcripts <- transcripts %>%
  mutate(across(all_of(expression_cols), ~ na_if(., "0"))) %>%
  mutate(across(all_of(expression_cols), ~ na_if(., "#DIV/0!"))) %>%
  mutate(across(all_of(expression_cols), as.numeric))

# Step 5: Calculate average knockdown ratio per transcript
transcripts <- transcripts %>%
  rowwise() %>%
  mutate(avg_knockdown = mean(c_across(all_of(expression_cols)), na.rm = TRUE)) %>%
  ungroup()

# Step 6: Bin average knockdown ratios
transcripts <- transcripts %>%
  mutate(knockdown_bin = case_when(
    avg_knockdown < 0.25 ~ "<0.25",
    avg_knockdown >= 0.25 & avg_knockdown <= 0.5 ~ "0.25–0.5",
    avg_knockdown > 0.5 & avg_knockdown <= 0.75 ~ "0.51–0.75",
    avg_knockdown > 0.75 & avg_knockdown <= 1.0 ~ "0.76–1.00",
    avg_knockdown > 1.0 ~ ">1.00",
    TRUE ~ NA_character_
  ))

# Step 7: Plot distribution of knockdown bins
ggplot(transcripts, aes(x = knockdown_bin)) +
  geom_bar(fill = "darkred", color = "black") +
  labs(
    title = "Distribution of Average Knockdown Ratios (Subset)",
    x = "Knockdown Ratio Bin",
    y = "Number of Transcripts"
  ) +
  theme_minimal()
