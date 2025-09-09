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
expr_df <- readRDS("Data/melanoma_transcript_cleaned.rds")
colnames(expr_df)[1] <- "transcript_id"

# Load transcript lists
spearman_ids <- read.csv("Data/TCGA_spearman_top1500.txt", stringsAsFactors = FALSE)
pearson_ids  <- read.csv("Data/TCGA_pearson_top_1500.txt", stringsAsFactors = FALSE)

# Combine and remove duplicates
combined_ids <- unique(rbind(spearman_ids, pearson_ids))

# Optional: ensure column is named 'transcript_id'
colnames(combined_ids) <- "transcript_id"

# Save combined list
write.csv(combined_ids, "Data/list_transcript.txt", row.names = FALSE)




# List of transcript IDs to keep
transcript_ids <- read.csv("Data/list_transcript.txt", stringsAsFactors = FALSE)

# Mapping file: transcript ID -> Gene
transcript_to_gene <- read.csv("Data/transcript_to_gene.csv", stringsAsFactors = FALSE)

# TCGA phenotype file
TCGA_phenotype <- read.csv("Data/TCGA_phenotype.csv", stringsAsFactors = FALSE)

# ---------------------------------------
# Step 2: Build Gene_transcript mapping
# ---------------------------------------

# Join transcript_ids with gene info
Gene_transcript <- inner_join(transcript_ids, transcript_to_gene, by = "transcript_id")

# -------------------------------------------
# Step 3: Subset Expression and Annotate Gene
# -------------------------------------------

# Subset expression to only selected transcripts
expr_subset <- inner_join(transcript_ids, expr_df, by = "transcript_id")
saveRDS(expr_subset, "expr_subset.rds")

# Merge gene names into expression data
expr_annotated <- left_join(expr_subset, Gene_transcript, by = "transcript_id") %>%
  select(transcript_id, Gene, everything())

# --------------------------------------------
# Step 4: Add primary_disease as a top row
# --------------------------------------------

sample_cols <- setdiff(colnames(expr_annotated), c("transcript_id", "Gene"))
phenotype_matched <- TCGA_phenotype[TCGA_phenotype$sample %in% sample_cols, ]
primary_disease_row <- phenotype_matched$primary_disease[match(sample_cols, phenotype_matched$sample)]
new_row <- c("primary_disease", NA, primary_disease_row)
new_row_df <- as.data.frame(t(new_row), stringsAsFactors = FALSE)
colnames(new_row_df) <- colnames(expr_annotated)
expr_annotated_with_disease <- rbind(new_row_df, expr_annotated)

# --------------------------------------
# Step 5: Reshape data for rank analysis
# --------------------------------------

expr_data <- expr_annotated_with_disease[-1, ]
expr_data <- expr_data %>% mutate(across(-c(transcript_id, Gene), as.numeric))

# Convert to long format
long_expr <- expr_data %>%
  pivot_longer(cols = -c(transcript_id, Gene), names_to = "sample", values_to = "expression") %>%
  left_join(TCGA_phenotype, by = "sample") %>%
  filter(!is.na(expression), !is.na(primary_disease))

# -------------------------------------------------------
# Step 6: Rank summary (mean + median) + Save table/plot
# -------------------------------------------------------

# Compute per-transcript ranks by disease
ranked_df <- long_expr %>%
  group_by(transcript_id, primary_disease) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(transcript_id) %>%
  mutate(rank = rank(-mean_expr, ties.method = "average")) %>%
  ungroup()

# Summarize mean/median rank
rank_summary <- ranked_df %>%
  group_by(primary_disease) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

# Save rank summary table
write.csv(rank_summary, "expression_rank_summary_table.csv", row.names = FALSE)
saveRDS(rank_summary, "expression_rank_summary_table.rds")

# Prepare for bar plot
rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value") %>%
  mutate(highlight = ifelse(primary_disease %in% c("SKCM", "UVM"), "highlight", "other"),
         fill_key = paste(statistic, highlight, sep = "."))

color_map <- c(
  "mean_rank.highlight" = "#ff1493",
  "median_rank.highlight" = "#8a2be2",
  "mean_rank.other" = "gray70",
  "median_rank.other" = "gray50"
)

# Plot rank summary and save to PDF
p_rank <- ggplot(rank_long, aes(x = reorder(primary_disease, rank_value), y = rank_value, fill = fill_key)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Mean and Median Rank of Expression by Disease",
    x = "Primary Disease",
    y = "Rank (lower = higher expression)",
    fill = "Statistic"
  ) +
  scale_fill_manual(
    values = color_map,
    labels = c(
      "mean_rank.highlight" = "Mean Rank (SKCM/UVM)",
      "mean_rank.other" = "Mean Rank (Other)",
      "median_rank.highlight" = "Median Rank (SKCM/UVM)",
      "median_rank.other" = "Median Rank (Other)"
    )
  )

ggsave("expression_rank_summary.pdf", plot = p_rank, width = 10, height = 6)
