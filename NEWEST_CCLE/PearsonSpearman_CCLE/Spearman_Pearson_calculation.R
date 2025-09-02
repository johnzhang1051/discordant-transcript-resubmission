library(data.table)
library(dplyr)
Gene <- fread("Cleaned_data_sets/Gene_expression_melanoma.csv")
Transcript <- fread("Cleaned_data_sets/Transcript_expression_melanoma_log2.csv")
### # Remove all rows where Sample_ID is "ACH-000931"as this is a duplicate in transcripts as has two different set of values
Transcript <- Transcript[Transcript$Sample_ID != "ACH-000931", ]
Gene <- Gene[Gene$Sample_ID != "ACH-000931", ]


transcript_type <-read.csv("Cleaned_data_sets/transcripttype.csv")


## MITF-M ENST00000394351
## rearrange data sets as needed to do Pearson and Spearman against MITF-M
# Required libraries
library(dplyr)
library(tibble)


# Assume your data frames are called Transcript and Gene
# Column 1 is Sample_ID, and the rest are expression values

# Ensure Sample_ID is character (not factor)
Transcript$Sample_ID <- as.character(Transcript[[1]])
Gene$Sample_ID <- as.character(Gene[[1]])

# Align by Sample_ID
common_samples <- intersect(Transcript$Sample_ID, Gene$Sample_ID)
Transcript_sub <- Transcript %>% filter(Sample_ID %in% common_samples)
Gene_sub <- Gene %>% filter(Sample_ID %in% common_samples)

# Drop Sample_ID for correlation and keep it separate
trans_expr <- Transcript_sub %>% column_to_rownames("Sample_ID")
gene_expr <- Gene_sub %>% column_to_rownames("Sample_ID")

# Check if reference transcript exists
ref_id <- "ENST00000394351"
if (!(ref_id %in% colnames(trans_expr))) {
  stop(paste("Reference transcript", ref_id, "not found in Transcript"))
}

# Extract reference expression vector
ref_vector <- trans_expr[[ref_id]]

# Function to compute correlations
get_cor_df <- function(mat, method = "spearman", ref_vector) {
  cor_vals <- apply(mat, 2, function(x) cor(ref_vector, x, method = method))
  data.frame(ID = colnames(mat), Correlation = cor_vals, row.names = NULL)
}

# Compute Spearman and Pearson for all transcripts
cor_trans_spearman <- get_cor_df(trans_expr, "spearman", ref_vector)
cor_trans_pearson  <- get_cor_df(trans_expr, "pearson",  ref_vector)


# View top correlations
head(cor_trans_spearman[order(-abs(cor_trans_spearman$Correlation)), ])


# Optional: Save to CSV
write.csv(cor_trans_spearman, "Transcript_Spearman_ENST00000394351.csv", row.names = FALSE)
write.csv(cor_trans_pearson,  "Transcript_Pearson_ENST00000394351.csv",  row.names = FALSE)
## subset protein coding and resave
# Ensure IDs are characters
transcript_type$ID <- as.character(transcript_type$ID)
cor_trans_spearman$ID <- as.character(cor_trans_spearman$ID)
cor_trans_pearson$ID  <- as.character(cor_trans_pearson$ID)

# Get only protein-coding transcript IDs
protein_coding_ids <- transcript_type %>%
  filter(transcript_type == "protein_coding") %>%
  pull(ID)

# Subset correlation lists
cor_trans_spearman_pc <- cor_trans_spearman %>%
  filter(ID %in% protein_coding_ids)

cor_trans_pearson_pc <- cor_trans_pearson %>%
  filter(ID %in% protein_coding_ids)

# Save filtered results
write.csv(cor_trans_spearman_pc, "cor_trans_spearman_protein_coding.csv", row.names = FALSE)
write.csv(cor_trans_pearson_pc,  "cor_trans_pearson_protein_coding.csv", row.names = FALSE)
# Remove rows with NA in Correlation column
cor_trans_spearman_pc_clean <- cor_trans_spearman_pc %>%
  filter(!is.na(Correlation))

cor_trans_pearson_pc_clean <- cor_trans_pearson_pc %>%
  filter(!is.na(Correlation))

# Save cleaned data frames
write.csv(cor_trans_spearman_pc_clean, "cor_trans_spearman_protein_coding_clean.csv", row.names = FALSE)
write.csv(cor_trans_pearson_pc_clean,  "cor_trans_pearson_protein_coding_clean.csv", row.names = FALSE)


##
# Ensure correct data types
Transcript$Sample_ID <- as.character(Transcript$Sample_ID)
Gene$Sample_ID <- as.character(Gene$Sample_ID)

# Intersect sample IDs
common_samples <- intersect(Transcript$Sample_ID, Gene$Sample_ID)

# Filter to matching samples only
Transcript_matched <- Transcript %>% filter(Sample_ID %in% common_samples)
Gene_matched <- Gene %>% filter(Sample_ID %in% common_samples)

# Confirm ordering is identical
Transcript_matched <- Transcript_matched %>% arrange(Sample_ID)
Gene_matched <- Gene_matched %>% arrange(Sample_ID)
# Convert to matrix with Sample_ID as rownames
trans_expr <- Transcript_matched %>% column_to_rownames("Sample_ID")
gene_expr <- Gene_matched %>% column_to_rownames("Sample_ID")

# Target transcript
target_id <- "ENST00000394351"
target_vector <- trans_expr[[target_id]]
# Initialize data frame
cor_results <- data.frame(
  Gene = colnames(gene_expr),
  Spearman = NA_real_,
  Pearson = NA_real_
)

# Loop through each gene
for (i in seq_along(cor_results$Gene)) {
  gene_vec <- gene_expr[[cor_results$Gene[i]]]
  
  # Check length match
  if (length(gene_vec) == length(target_vector)) {
    cor_results$Spearman[i] <- cor(target_vector, gene_vec, method = "spearman")
    cor_results$Pearson[i]  <- cor(target_vector, gene_vec, method = "pearson")
  }
}
# Sort by descending correlation
cor_results_sorted <- cor_results %>%
  arrange(desc(Spearman))

# Optional: View top results
head(cor_results_sorted, 10)



