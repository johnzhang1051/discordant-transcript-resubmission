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

library(dplyr)
library(readr)

# expr_data: assumes already loaded, with columns: transcript_id, Gene, sample_1, sample_2, ...

# Get sample-to-disease mapping
sample_map <- TCGA_phenotype %>% select(sample, primary_disease)

# Keep only the samples present in the expression data
valid_samples <- intersect(colnames(expr_data), sample_map$sample)

# Subset expression matrix to samples
expr_matrix <- expr_data[, c("transcript_id", "Gene", valid_samples)]

# Transpose for long format manually, then group by disease
expr_long_list <- lapply(valid_samples, function(samp) {
  data.frame(
    transcript_id = expr_matrix$transcript_id,
    sample = samp,
    expression = expr_matrix[[samp]],
    stringsAsFactors = FALSE
  )
})

# Bind rows (still large, but controlled)
expr_long_df <- bind_rows(expr_long_list)

# Merge with disease labels
expr_long_df <- expr_long_df %>%
  left_join(sample_map, by = "sample") %>%
  filter(!is.na(expression), !is.na(primary_disease))

# Now compute rank summary
ranked_df <- expr_long_df %>%
  group_by(transcript_id, primary_disease) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(transcript_id) %>%
  mutate(rank = rank(-mean_expr, ties.method = "average")) %>%
  ungroup()

rank_summary <- ranked_df %>%
  group_by(primary_disease) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

# Save summary
write.csv(rank_summary, "expression_rank_summary_table.csv", row.names = FALSE)
saveRDS(rank_summary, "expression_rank_summary_table.rds")

# -----------------------------
# Plot mean/median rank summary
# -----------------------------

library(ggplot2)
library(forcats)

# Reshape to long format for plotting
rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value") %>%
  mutate(
    highlight = ifelse(primary_disease %in% c("SKCM", "UVM"), "highlight", "other"),
    fill_key = paste(statistic, highlight, sep = ".")
  )

# Define fill colors
color_map <- c(
  "mean_rank.highlight" = "#ff1493",    # SKCM/UVM mean rank - bright pink
  "median_rank.highlight" = "#8a2be2",  # SKCM/UVM median rank - bright purple
  "mean_rank.other" = "gray70",
  "median_rank.other" = "gray50"
)

# Create bar plot
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

# Save plot to PDF
ggsave("expression_rank_summary.pdf", plot = p_rank, width = 10, height = 6)

