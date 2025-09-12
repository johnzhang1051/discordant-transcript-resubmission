# Load libraries
library(data.table)
library(dplyr)

# Load raw data
expr_df <- read.csv("Data/CCLE_expression.csv")
model_ID <- read.csv("data/cell_line_disease.csv")

# Rename the first column
colnames(expr_df)[1] <- "ModelID"

# Merge with model_ID data frame by "ModelID"
merged_df <- merge(model_ID, expr_df, by = "ModelID")

# Clean column names (remove the ".." patterns)
colnames(merged_df)[4:ncol(merged_df)] <- gsub("\\.\\..*?\\.", "", colnames(merged_df)[4:ncol(merged_df)])

# Add gene filtering
gene_list <- read.csv("data/final_gene_list.txt")
filtered_df <- merged_df[, c("ModelID", "CellLineName", "DepmapModelType", intersect(gene_list$Gene, colnames(merged_df)))]

# TCGA DepMap Filtering to "allowed" DepMap Model Types
TCGA_model_type <- read.csv("data/TCGA_DepmapModelType.txt")
TCGA_modeltype_unique <- unique(trimws(as.character(TCGA_model_type$DepmapModelType)))
filtered_df$DepmapModelType <- trimws(as.character(filtered_df$DepmapModelType))

# Filter to keep only matching model types
filtered_df <- filtered_df[filtered_df$DepmapModelType %in% TCGA_modeltype_unique, ]
# Save with both gene and model type filtering
saveRDS(filtered_df, "data/CCLE_gene_disease_Figure1B_new.rds")
