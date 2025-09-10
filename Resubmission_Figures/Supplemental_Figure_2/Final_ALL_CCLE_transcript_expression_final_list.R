# Load libraries
library(data.table)
library(dplyr)

# Load raw data
expr_df <- readRDS("data/CCLE_transcript_expression.rds")
model_ID <- read.csv("data/cell_line_disease.csv")

# Extract column names (skip first column)
col_names <- colnames(expr_df)[-1]

# Split on double period to get gene and transcript info
split_names <- strsplit(col_names, "\\.\\.")
gene_names <- sapply(split_names, `[`, 1)
transcript_ids <- sapply(split_names, `[`, 2)

# Create header rows
header1 <- c(NA, gene_names)
header2 <- c(NA, transcript_ids)

# Combine all as a data frame
full_df <- rbind(header1, header2, as.matrix(expr_df))
full_df <- as.data.frame(full_df, stringsAsFactors = FALSE)

# Set transcript IDs as column names
colnames(full_df) <- transcript_ids
colnames(full_df)[1] <- "ModelID"

# Remove the first two rows (headers)
expr_temp_data <- full_df[-c(1, 2), ]
expr_temp_data$ModelID <- as.character(expr_temp_data$ModelID)

# Convert expression columns to numeric
expr_temp_data[, -1] <- lapply(expr_temp_data[, -1], as.numeric)

# Merge with model metadata
merged_expr <- merge(model_ID, expr_temp_data, by = "ModelID")

# Clean column names
colnames(merged_expr)[4:ncol(merged_expr)] <- sub("\\.$", "", colnames(merged_expr)[4:ncol(merged_expr)])

# Add transcript filtering
final_list <- read.csv("data/final_transcript_list.txt")
final_transcripts <- as.character(final_list$transcript_id)
final_transcripts <- sub("\\.$", "", trimws(final_transcripts))

# Clean column names for matching
colnames(merged_expr) <- sub("\\.$", "", trimws(colnames(merged_expr)))

# Find matching transcript columns
matching_cols <- which(colnames(merged_expr) %in% final_transcripts)

# Keep metadata (cols 1-3) and matching expression columns
filtered_df <- merged_expr[, c(1:3, matching_cols[matching_cols > 3])]

# Save with both model type and transcript filtering
saveRDS(filtered_df, "data/Figure_1C_CCLE_transcript_new.rds")
