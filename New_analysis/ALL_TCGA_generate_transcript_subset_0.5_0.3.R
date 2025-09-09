library(data.table)

expr_df <- readRDS("Data/melanoma_transcript_cleaned.rds")
expr_df[1:5, 1:5]

# Step 1: Read your list of transcript IDs
transcript_ids <- read.csv("Data/correlation_only_list.csv")


# Step 3: Rename the first column to 'transcript_id'
colnames(expr_df)[1] <- "transcript_id"

# Step 4: Perform an inner join to keep only matching transcript_ids
library(dplyr)

expr_subset <- inner_join(transcript_ids, expr_df, by = "transcript_id")

# Step 5: Save the subset if desired
saveRDS(expr_subset, "expr_subset.rds")

TCGA_phenotype <- read.csv("Data/TCGA_phenotype.csv")
Gene_transcript <- read.csv("Data/final_transcript_list.csv")

# Assuming both data frames are already loaded
# expr_subset: contains 'transcript_id' + expression values
# Gene_transcript: has 'transcript_id' and 'gene_name'

library(dplyr)

# Merge to add gene_name column to expr_subset
expr_annotated <- left_join(expr_subset, Gene_transcript, by = "transcript_id")

expr_annotated <- expr_annotated %>%
  select(transcript_id, Gene, everything())

# Step 1: Identify which columns in expr_annotated are sample IDs
sample_cols <- setdiff(colnames(expr_annotated), c("transcript_id", "Gene"))

# Step 2: Filter phenotype data to just those sample IDs
phenotype_matched <- TCGA_phenotype[TCGA_phenotype$sample %in% sample_cols, ]

# Step 3: Create a named vector of primary diseases for those samples
# Match in the correct order of sample columns
primary_disease_row <- phenotype_matched$primary_disease[match(sample_cols, phenotype_matched$sample)]

# Step 4: Build the new first row (with labels for col 1 and col 2)
new_row <- c("primary_disease", NA, primary_disease_row)

# Step 5: Turn it into a one-row data frame with the same column names
new_row_df <- as.data.frame(t(new_row), stringsAsFactors = FALSE)
colnames(new_row_df) <- colnames(expr_annotated)

# Step 6: Bind the new row to the top of expr_annotated
expr_annotated_with_disease <- rbind(new_row_df, expr_annotated)


# Remove the disease row
expr_data <- expr_annotated_with_disease[-1, ]

# Convert expression columns to numeric
expr_data <- expr_data %>%
  mutate(across(-c(transcript_id, Gene), as.numeric))

# Pivot to long format
long_expr <- expr_data %>%
  pivot_longer(
    cols = -c(transcript_id, Gene),
    names_to = "sample",
    values_to = "expression"
  )

# Join with phenotype
phenotype_subset <- TCGA_phenotype %>% select(sample, primary_disease)
long_expr <- left_join(long_expr, phenotype_subset, by = "sample")

# Remove missing values
long_expr <- long_expr %>%
  filter(!is.na(expression), !is.na(primary_disease))

# Color palette
unique_diseases <- unique(long_expr$primary_disease)
disease_palette <- setNames(rep("gray70", length(unique_diseases)), unique_diseases)
disease_palette["SKCM"] <- "#ff1493"  # Bright pink
disease_palette["UVM"] <- "#8a2be2"   # Bright purple
other_diseases <- setdiff(unique_diseases, c("SKCM", "UVM"))
safe_colors <- RColorBrewer::brewer.pal(min(length(other_diseases), 8), "Set2")
disease_palette[other_diseases] <- safe_colors

# Multi-page PDF output
pdf("transcript_expression_by_disease.pdf", width = 10, height = 6)

transcript_ids <- unique(long_expr$transcript_id)

for (tid in transcript_ids) {
  
  df <- long_expr %>% filter(transcript_id == tid)
  if (nrow(df) == 0) next
  
  df <- df %>%
    group_by(primary_disease) %>%
    mutate(median_expr = median(expression, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(primary_disease = fct_reorder(primary_disease, median_expr))
  
  gene_name <- unique(df$Gene)
  
  p <- ggplot(df, aes(x = primary_disease, y = expression, fill = primary_disease)) +
    geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5) +
    scale_fill_manual(values = disease_palette) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = "Primary Disease",
      y = "Expression (TPM)",
      fill = "Disease Type",
      title = paste0("Gene: ", gene_name, "\nTranscript ID: ", tid)
    )
  
  print(p)
}

dev.off()

## plot ranks

library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Compute mean expression per transcript × disease
ranked_df <- long_expr %>%
  group_by(transcript_id, primary_disease) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(transcript_id) %>%
  mutate(rank = rank(-mean_expr, ties.method = "average")) %>%  # Rank 1 = highest expression
  ungroup()

# Step 2: Calculate both mean and median rank per disease
rank_summary <- ranked_df %>%
  group_by(primary_disease) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Pivot longer for ggplot
rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value")

# Step 4: Plot
ggplot(rank_long, aes(x = reorder(primary_disease, rank_value), y = rank_value, fill = statistic)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Mean and Median Rank of Expression by Disease",
    x = "Primary Disease",
    y = "Rank (lower = higher expression)",
    fill = "Statistic"
  ) +
  scale_fill_manual(values = c("mean_rank" = "#3182bd", "median_rank" = "#e6550d"))

## color pallette
library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Compute mean expression per transcript × disease
ranked_df <- long_expr %>%
  group_by(transcript_id, primary_disease) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(transcript_id) %>%
  mutate(rank = rank(-mean_expr, ties.method = "average")) %>%  # Rank 1 = highest expression
  ungroup()

# Step 2: Calculate both mean and median rank per disease
rank_summary <- ranked_df %>%
  group_by(primary_disease) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Pivot longer for plotting
rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value")

# Step 4: Add a highlight variable
rank_long <- rank_long %>%
  mutate(highlight = ifelse(primary_disease %in% c("SKCM", "UVM"), "highlight", "other"))

# Step 5: Define custom colors
color_map <- c(
  "mean_rank.highlight" = "#ff1493",  # bright pink for SKCM/UVM (mean)
  "median_rank.highlight" = "#8a2be2", # bright purple for SKCM/UVM (median)
  "mean_rank.other" = "gray70",
  "median_rank.other" = "gray50"
)

# Step 6: Build the fill key
rank_long <- rank_long %>%
  mutate(fill_key = paste(statistic, highlight, sep = "."))

# Step 7: Plot and save
p <- ggplot(rank_long, aes(x = reorder(primary_disease, rank_value), y = rank_value, fill = fill_key)) +
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

# Save to PDF
ggsave("expression_rank_summary.pdf", plot = p, width = 10, height = 6)




