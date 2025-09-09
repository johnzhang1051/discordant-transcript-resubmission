# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(cowplot)

# === 1. LOAD DATA ===
expr_df <- readRDS("Data/CCLE_transcript_expression.rds")
# Extract column names (skip first column, 'X')
col_names <- colnames(expr_df)[-1]

# Split on double period
split_names <- strsplit(col_names, "\\.\\.")

# Extract gene and transcript
gene_names <- sapply(split_names, `[`, 1)
transcript_ids <- sapply(split_names, `[`, 2)

# Create header rows
header1 <- c(NA, gene_names)
header2 <- c(NA, transcript_ids)

# Grab the data body (including sample names in first column)
data_body <- expr_df

# Combine all as a data frame
full_df <- rbind(header1, header2, as.matrix(data_body))
full_df <- as.data.frame(full_df, stringsAsFactors = FALSE)

# Optional: set meaningful column names (e.g., V1, V2, ...)
colnames(full_df) <- c("Sample", paste0("V", 2:ncol(full_df)))

# Save as RDS
saveRDS(full_df, file = "CCLE_transcript_needs_cleaned.rds")

### START HERE
expr_temp <- readRDS("Data/CCLE_transcript_needs_cleaned.rds")
colnames(expr_temp)[1] <- "ModelID"
model_ID <-  read.csv("Data/cell_line_disease.csv")

# Step 1: Extract gene names and transcript IDs
gene_names <- as.character(expr_temp[1, ])
transcript_ids <- as.character(expr_temp[2, ])

# Step 2: Set transcript IDs as column names
colnames(expr_temp) <- transcript_ids
colnames(expr_temp)[1] <- "ModelID"  # label first column properly

# Step 3: Remove the first two rows (we've saved them)
expr_temp_data <- expr_temp[-c(1, 2), ]

# Step 4: Convert to data frame and fix ModelID column
expr_temp_data <- as.data.frame(expr_temp_data)
expr_temp_data$ModelID <- as.character(expr_temp_data$ModelID)
colnames(expr_temp_data) <- transcript_ids
colnames(expr_temp_data)[1] <- "ModelID"

# Step 5: Convert expression columns to numeric
expr_temp_data[ , -1] <- lapply(expr_temp_data[ , -1], as.numeric)

# Step 6: Merge with model_id
merged_expr <- merge(model_ID, expr_temp_data, by = "ModelID")
colnames(merged_expr)[4:198840] <- sub("\\.$", "", colnames(merged_expr)[4:198840])


# Step 7: Create transcript → gene map
transcript_gene_map <- data.frame(
  TranscriptID = transcript_ids[-1],  # exclude "ModelID"
  Gene = gene_names[-1],
  stringsAsFactors = FALSE
)






TCGA_modeltype <- read.csv("Data/TCGA_DepmapModelType.txt")
TCGA_modeltype_unique <- unique(trimws(as.character(TCGA_modeltype$DepmapModelType)))
# Clean up the column for matching
merged_expr$DepmapModelType <- trimws(as.character(merged_expr$DepmapModelType))
# Filter to keep only matching model types
filtered_CCLE_transcript <- merged_expr[merged_expr$DepmapModelType %in% TCGA_modeltype_unique, ]

### multiple lists


# Define your input list files and short labels
transcript_lists <- list(
  list1 = "Data/final_transcript_list.txt",
  list2 = "Data/correlation_only_list.txt"
)

# Create storage for summaries
combined_summary <- list()
combined_ranks <- list()

for (list_name in names(transcript_lists)) {
  message("Processing ", list_name, "...")
  
  # === Load transcript list ===
  final_list <- read.csv(transcript_lists[[list_name]])
  final_transcripts <- sub("\\.$", "", trimws(as.character(final_list$transcript_id)))
  
  # Clean expression column names
  colnames(filtered_CCLE_transcript) <- sub("\\.$", "", trimws(colnames(filtered_CCLE_transcript)))
  
  # Subset filtered expression data
  matching_cols <- which(colnames(filtered_CCLE_transcript) %in% final_transcripts)
  filtered_df <- filtered_CCLE_transcript[, c(1:3, matching_cols[matching_cols > 3])]
  
  # Gather into long format
  long_df <- filtered_df %>%
    pivot_longer(cols = 4:ncol(filtered_df), names_to = "Transcript", values_to = "Expression")
  
  # Per-transcript rank stats
  rank_stats_df <- long_df %>%
    group_by(DepmapModelType, Transcript) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE),
              median_expr = median(Expression, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Color = ifelse(DepmapModelType == "SKCM", "pink", "gray"),
           GeneSet = list_name)
  
  # Mean + Median summary per model type
  summary_means <- rank_stats_df %>%
    group_by(DepmapModelType, Color, GeneSet) %>%
    summarise(Value = mean(mean_expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Mean")
  
  summary_medians <- rank_stats_df %>%
    group_by(DepmapModelType, Color, GeneSet) %>%
    summarise(Value = mean(median_expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Median")
  
  combined_summary[[list_name]] <- bind_rows(summary_means, summary_medians)
  
  # Rank-based summarization
  ranked_df <- rank_stats_df %>%
    group_by(Transcript) %>%
    mutate(Rank = rank(-mean_expr, ties.method = "average"))  # Rank across model types
  
  average_ranks <- ranked_df %>%
    group_by(DepmapModelType, GeneSet) %>%
    summarise(AverageRank = mean(Rank), .groups = "drop") %>%
    mutate(Color = ifelse(DepmapModelType == "SKCM", "pink", "gray"))
  
  combined_ranks[[list_name]] <- average_ranks
  
  # Optionally save intermediate files:
  write.csv(rank_stats_df, paste0("RankSummaries/RankStats_", list_name, ".csv"), row.names = FALSE)
}

# Combine for cross-list comparison
summary_all <- bind_rows(combined_summary)
ranks_all <- bind_rows(combined_ranks)

# === Combined Mean & Median Plot ===
# Combined Mean & Median Plot
combined_plot <- ggplot(summary_all, aes(x = DepmapModelType, y = Value, fill = Color)) +
  geom_col(position = position_dodge(width = 0.8), aes(group = Metric)) +
  facet_grid(Metric ~ GeneSet, scales = "free_y") +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean & Median Transcript Expression Across Gene Sets",
       x = "DepMap Model Type", y = "Expression")

ggsave("RankSummaries/Expression_Comparison_AllGeneSets.pdf", combined_plot, width = 12, height = 6)


# Rank-Based Summary Plot
rank_plot <- ggplot(ranks_all, aes(x = fct_reorder(DepmapModelType, AverageRank), y = AverageRank, fill = Color)) +
  geom_col() +
  facet_wrap(~ GeneSet) +
  scale_fill_identity() +
  theme_bw() +
  labs(title = "Average Rank of Transcript Expression Across Gene Sets",
       x = "DepMap Model Type", y = "Average Rank (Lower = Higher)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("RankSummaries/AverageRank_Comparison_AllGeneSets.pdf", rank_plot, width = 10, height = 6)


# Compare SKCM expression across gene sets
skcm_expr <- summary_all %>%
  filter(DepmapModelType == "SKCM", Metric == "Mean")  # or "Median"

# ANOVA or Kruskal-Wallis across gene sets
kruskal.test(Value ~ GeneSet, data = skcm_expr)


### ENRICHMENT
# Assume you already have: long_df with Transcript, DepmapModelType, Expression
# Create storage for enrichment comparison
enrichment_all <- list()

for (list_name in names(transcript_lists)) {
  message("Processing enrichment for ", list_name, "...")
  
  # Load list and clean
  final_list <- read.csv(transcript_lists[[list_name]])
  final_transcripts <- sub("\\.$", "", trimws(as.character(final_list$transcript_id)))
  
  # Filter expression data
  matching_cols <- which(colnames(filtered_CCLE_transcript) %in% final_transcripts)
  filtered_df <- filtered_CCLE_transcript[, c(1:3, matching_cols[matching_cols > 3])]
  
  # Long format
  long_df <- filtered_df %>%
    pivot_longer(cols = 4:ncol(filtered_df), names_to = "Transcript", values_to = "Expression")
  
  # Compute enrichment for SKCM vs Others
  enrichment_df <- long_df %>%
    group_by(Transcript, DepmapModelType) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = ifelse(DepmapModelType == "SKCM", "SKCM", "Other")) %>%
    group_by(Transcript, Group) %>%
    summarise(group_mean = mean(mean_expr), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = group_mean) %>%
    mutate(log2FC_SKCM = log2(SKCM / Other),
           GeneSet = list_name)
  
  enrichment_all[[list_name]] <- enrichment_df
}

# Combine all
enrichment_combined <- bind_rows(enrichment_all)

# Plot
ggplot(enrichment_combined, aes(x = GeneSet, y = log2FC_SKCM, fill = GeneSet)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "SKCM Enrichment (log2 FC of SKCM vs Other)",
       y = "log2 Fold-Change (SKCM / Other)",
       x = "Transcript List")

# Statistical test
ttest_result <- t.test(log2FC_SKCM ~ GeneSet, data = enrichment_combined)
print(ttest_result)



