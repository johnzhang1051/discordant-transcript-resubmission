# Load required libraries
library(tidyverse)

# Step 1: Load the annotated transcripts
transcripts <- readRDS("Data/hg38_annotated_transcripts_for_analysis.rds")

# Make sure transcript_id is a column (not row names)
if (is.null(transcripts$transcript_id)) {
  transcripts <- transcripts %>%
    rownames_to_column(var = "transcript_id")
}



# Ensure MITF_peak_n and EBOX_n are numeric
transcripts <- transcripts %>%
  mutate(
    MITF_peak_n = as.numeric(MITF_peak_n),
    EBOX_n = as.numeric(EBOX_n)
  )


# Step 2: Create and display plots



# Plot 1: MITF_peak_n
ggplot(transcripts, aes(x = MITF_peak_n)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(
    title = "Distribution of MITF Peaks per Transcript",
    x = "MITF_peak_n",
    y = "Count"
  ) +
  theme_minimal()

# Plot 2: EBOX_n
ggplot(transcripts, aes(x = EBOX_n)) +
  geom_histogram(binwidth = 1, fill = "darkgreen", color = "black") +
  scale_x_continuous(breaks = 0:9, limits = c(0, 9)) +
  labs(
    title = "Distribution of EBOX Count per Transcript",
    x = "EBOX_n",
    y = "Count"
  ) +
  theme_minimal()


# Plot 3: Remote_TSS (all transcripts)
ggplot(transcripts, aes(x = Remote_TSS)) +
  geom_bar(fill = "purple", color = "black") +
  labs(
    title = "Remote TSS Classification (All Transcripts)",
    x = "Remote_TSS",
    y = "Count"
  ) +
  theme_minimal()

# Plot 4: Remote_TSS (protein_coding only)
transcripts_protein_coding <- transcripts %>%
  filter(transcript_type == "protein_coding")

ggplot(transcripts_protein_coding, aes(x = Remote_TSS)) +
  geom_bar(fill = "orange", color = "black") +
  labs(
    title = "Remote TSS Classification (Protein-Coding Transcripts)",
    x = "Remote_TSS",
    y = "Count"
  ) +
  theme_minimal()


## viewing some transcripts
# Load required libraries
library(tidyverse)

# Step 1: Import the list of transcript IDs
final_transcripts <- read_tsv("Data/final_transcript_list.txt", col_names = "transcript_id")

# Step 2: Filter the transcripts object based on all three conditions
selected_transcripts <- transcripts %>%
  filter(
    transcript_id %in% final_transcripts$transcript_id,
    transcript_type == "protein_coding",
    Remote_TSS == FALSE
  )

# Step 3: View or inspect
glimpse(selected_transcripts)
# Or just the first few:
head(selected_transcripts)


