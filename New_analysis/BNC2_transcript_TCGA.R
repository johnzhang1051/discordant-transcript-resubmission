# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(readr)

# -------------------------
# Step 1: Load Input Files
# -------------------------

# Expression matrix
expr_df <- readRDS("Data/TCGA_transcript_cleaned.rds")
colnames(expr_df)[1] <- "transcript_id"


# List of transcript IDs to keep ***NEEDS TO BE UPDATED WITH A CUSTOM LIST***
transcript_ids <- read.csv("Data/BNC2_list_transcript.txt", stringsAsFactors = FALSE)
# TCGA phenotype file
TCGA_phenotype <- read.csv("Data/TCGA_phenotype.csv", stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Filter expr_df to only include selected transcript_ids
expr_filtered <- expr_df %>%
  filter(transcript_id %in% transcript_ids$transcript_id)

# Step 2: Convert to long format (transcript_id, sample, expression)
expr_long <- expr_filtered %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "sample",
    values_to = "expression"
  )

# Step 3: Merge with phenotype
expr_annotated <- expr_long %>%
  left_join(TCGA_phenotype, by = "sample") %>%
  filter(!is.na(primary_disease))

# Step 4: Loop through each transcript and plot separately
unique_transcripts <- unique(expr_annotated$transcript_id)

for (tx in unique_transcripts) {
  plot_data <- expr_annotated %>% filter(transcript_id == tx)
  
  p <- ggplot(plot_data, aes(x = primary_disease, y = expression, fill = primary_disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    coord_cartesian(ylim = c(-5, 10)) +  # You can modify limits if needed
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste("Expression of", tx, "by Primary Disease"),
      x = "Primary Disease",
      y = "Expression"
    )
  
  ggsave(filename = paste0("expression_", tx, ".pdf"), plot = p, width = 8, height = 6)
}
