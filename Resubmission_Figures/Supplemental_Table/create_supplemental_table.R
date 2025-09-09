# Minimal imports - dplyr includes most of what we need
library(dplyr)
library(data.table)
library(tidyr)

# Load all necessary data files:

# CCLE correlation results (both Pearson and Spearman)
ccle_pearson <- read.csv("data/CoCor_analysis/cocor_all_results.csv")
ccle_spearman <- read.csv("data/CoCor_analysis/cocor_spearman_all_results.csv")

# Tsoi correlation results
tsoi_correlations <- read.csv("data/Tsoi/correlation_to_ENST00000394351.csv")

# Transcript type information
transcript_type <- read.csv("data/transcripttype.csv")

# Discordant and correlated transcript lists
discordant_list <- read.csv("data/final_paper_lists/discordant_RESUBMISSION.csv")
correlated_list <- read.csv("data/final_paper_lists/correlated_RESUBMISSION.csv")

# Transcript expression for filtering
transcript_expression <- fread("data/Transcript_expression_melanoma.csv")

# Additional data sources (add paths as needed)
protein_coding_filtered <- fread("data/final_paper_lists/protein_coding_RESUBMISSION.csv")

# chip_peaks <- read.csv("path/to/chip_seq_peaks.csv")  # if available
# transcript_coords <- read.csv("path/to/transcript_coordinates.csv")  # if available

# Step 1: Start with all transcripts from correlation analyses
# Combine CCLE results (both Pearson and Spearman have same transcript coverage)
all_transcripts <- unique(c(ccle_pearson$transcript_ID, tsoi_correlations$transcript_id))

# Create base dataset
comprehensive_dataset <- data.frame(
  transcript_id = all_transcripts,
  stringsAsFactors = FALSE
)

# Step 2: Add CCLE correlation data
# Merge Pearson results
ccle_pearson_clean <- ccle_pearson %>%
  select(transcript_id = transcript_ID, 
         Gene,
         mitf_pearson_correlation = r_transcript,
         mitf_gene_pearson_correlation = r_gene,
         ccle_pearson_fdr = FDR,
         ccle_pearson_pvalue = p_value)

comprehensive_dataset <- left_join(comprehensive_dataset, ccle_pearson_clean, by = "transcript_id")

# Merge Spearman results
ccle_spearman_clean <- ccle_spearman %>%
  select(transcript_id = transcript_ID,
         mitf_spearman_correlation = r_transcript,
         mitf_gene_spearman_correlation = r_gene,
         ccle_spearman_fdr = FDR,
         ccle_spearman_pvalue = p_value)

comprehensive_dataset <- left_join(comprehensive_dataset, ccle_spearman_clean, by = "transcript_id")

# Step 3: Add Tsoi correlation data
tsoi_clean <- tsoi_correlations %>%
  select(transcript_id,
         tsoi_pearson_correlation = Pearson,
         tsoi_spearman_correlation = Spearman)

comprehensive_dataset <- left_join(comprehensive_dataset, tsoi_clean, by = "transcript_id")

# Step 4: Add transcript type information
transcript_type_clean <- transcript_type %>%
  select(transcript_id = transcript_ID, transcript_type)

comprehensive_dataset <- left_join(comprehensive_dataset, transcript_type_clean, by = "transcript_id")

# Step 5: Add gene_id (extract from transcript_id or use mapping)
# If you have a transcript-to-gene mapping file, use that
# Otherwise, extract gene_id from existing Gene column or transcript_id
comprehensive_dataset$gene_id <- comprehensive_dataset$Gene  # Use Gene column from CCLE data

# Step 6: Add data source information
comprehensive_dataset$data_source <- case_when(
  !is.na(comprehensive_dataset$ccle_pearson_fdr) & !is.na(comprehensive_dataset$tsoi_pearson_correlation) ~ "CCLE+Tsoi",
  !is.na(comprehensive_dataset$ccle_pearson_fdr) ~ "CCLE",
  !is.na(comprehensive_dataset$tsoi_pearson_correlation) ~ "Tsoi",
  TRUE ~ "Unknown"
)


# Step 9: Add classification flags
comprehensive_dataset$is_discordant <- comprehensive_dataset$transcript_id %in% discordant_list$transcript_id
comprehensive_dataset$is_correlated <- comprehensive_dataset$transcript_id %in% correlated_list$transcript_id
comprehensive_dataset$is_protein_coding_filtered <- comprehensive_dataset$transcript_id %in% protein_coding_filtered$transcript_id

# Step 10: Add combined FDR (use the more significant one from CCLE analyses)
comprehensive_dataset$fdr <- pmin(comprehensive_dataset$ccle_pearson_fdr, 
                                  comprehensive_dataset$ccle_spearman_fdr, 
                                  na.rm = TRUE)

# Step 11: Add count information (mean expression across samples)
count_info <- transcript_expression_clean %>%
  summarise_at(vars(-Sample_ID), mean, na.rm = TRUE) %>%
  pivot_longer(everything(), names_to = "transcript_id", values_to = "mean_counts")

comprehensive_dataset <- left_join(comprehensive_dataset, count_info, by = "transcript_id")
comprehensive_dataset$counts <- comprehensive_dataset$mean_counts

# Step 12: Add ChIP-seq peak information
# Load and process ChIP-seq data (from your existing code)
Kenny <- read.csv("data/chip/Kenny.csv")
Laurette <- read.csv("data/chip/Laurette.csv")
Louph <- read.csv("data/chip/Louphrasitthiphol.csv")

# Standardize transcript IDs
standardize_transcripts <- function(df) {
  df %>%
    dplyr::rename(transcript_id = transcriptId) %>%
    dplyr::mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))
}

Kenny <- standardize_transcripts(Kenny)
Laurette <- standardize_transcripts(Laurette)
Louph <- standardize_transcripts(Louph)

# Filter for promoter peaks (<=1kb)
filter_promoters <- function(df) {
  df %>%
    dplyr::filter(annotation == "Promoter (<=1kb)") %>%
    dplyr::select(transcript_id, annotation)
}

Kenny_promoters <- filter_promoters(Kenny)
Laurette_promoters <- filter_promoters(Laurette)
Louph_promoters <- filter_promoters(Louph)

# Combine all peaks (any dataset)
all_chip_peaks <- unique(c(
  Kenny_promoters$transcript_id,
  Laurette_promoters$transcript_id,
  Louph_promoters$transcript_id
))

# Add ChIP-seq peak information
comprehensive_dataset$has_chip_peak <- comprehensive_dataset$transcript_id %in% all_chip_peaks

# Add specific peak counts per dataset
comprehensive_dataset$kenny_peak <- comprehensive_dataset$transcript_id %in% Kenny_promoters$transcript_id
comprehensive_dataset$laurette_peak <- comprehensive_dataset$transcript_id %in% Laurette_promoters$transcript_id
comprehensive_dataset$louph_peak <- comprehensive_dataset$transcript_id %in% Louph_promoters$transcript_id
comprehensive_dataset$total_chip_peaks <- as.integer(comprehensive_dataset$kenny_peak) + 
  as.integer(comprehensive_dataset$laurette_peak) + 
  as.integer(comprehensive_dataset$louph_peak)

# Step 13: Add MITF Overexpression values
# Load overexpression datasets
GSE163646 <- read.csv("data/mitf_oe/GSE_163646_OE_kallisto_transcript_TPM.csv")
PRJNA704810 <- read.csv("data/mitf_oe/PRJNA704810_OE_kallisto_transcript_TPM.csv")

# Merge datasets
merged_OE <- merge(GSE163646, PRJNA704810, by="transcript_id")

# Define sample groups
prjna_con <- c("SRR13782518", "SRR13782519")
prjna_oe  <- c("SRR13782520", "SRR13782521", "SRR13782522")
gse_oe    <- c("SRR13282354", "SRR13282355", "SRR13282356")
gse_con   <- c("SRR13282357", "SRR13282358", "SRR13282359")

# Calculate group means
mean_expr_df <- merged_OE %>%
  mutate(
    PRJNA_CON_mean = rowMeans(select(., all_of(prjna_con)), na.rm = TRUE),
    PRJNA_OE_mean  = rowMeans(select(., all_of(prjna_oe)),  na.rm = TRUE),
    GSE_CON_mean   = rowMeans(select(., all_of(gse_con)),   na.rm = TRUE),
    GSE_OE_mean    = rowMeans(select(., all_of(gse_oe)),    na.rm = TRUE)
  ) %>%
  select(transcript_id, PRJNA_CON_mean, PRJNA_OE_mean, GSE_CON_mean, GSE_OE_mean)

# Calculate OE/CON ratios
ratio_df <- mean_expr_df %>%
  mutate(
    PRJNA_OE_CON_ratio = PRJNA_OE_mean / PRJNA_CON_mean,
    GSE_OE_CON_ratio   = GSE_OE_mean / GSE_CON_mean
  ) %>%
  select(transcript_id, PRJNA_OE_CON_ratio, GSE_OE_CON_ratio) %>%
  filter(is.finite(PRJNA_OE_CON_ratio), is.finite(GSE_OE_CON_ratio))

# Create simplified overexpression metric (average of both datasets)
overexpression_data <- ratio_df %>%
  mutate(
    mitf_overexpression = (PRJNA_OE_CON_ratio + GSE_OE_CON_ratio) / 2,
    prjna_overexpression = PRJNA_OE_CON_ratio,
    gse_overexpression = GSE_OE_CON_ratio
  ) %>%
  select(transcript_id, mitf_overexpression, prjna_overexpression, gse_overexpression)

# Merge with comprehensive dataset
comprehensive_dataset <- left_join(comprehensive_dataset, overexpression_data, by = "transcript_id")
# Step 13: Clean up and reorder columns
final_columns <- c(
  "transcript_id",
  "data_source",
  "gene_id",
  "transcript_type",
  "mitf_pearson_correlation",
  "mitf_spearman_correlation", 
  "mitf_gene_pearson_correlation",
  "mitf_gene_spearman_correlation",
  "fdr",
  "counts",
  "prjna_overexpression",
  "gse_overexpression",  "has_chip_peak",
  "is_correlated",
  "is_discordant",
  "is_protein_coding_filtered"
)

# Select and reorder columns
comprehensive_dataset_final <- comprehensive_dataset %>%
  select(all_of(final_columns[final_columns %in% colnames(comprehensive_dataset)]))

# Add any missing columns as NA
missing_cols <- setdiff(final_columns, colnames(comprehensive_dataset_final))
for (col in missing_cols) {
  comprehensive_dataset_final[[col]] <- NA
}

# Reorder to match desired order
comprehensive_dataset_final <- comprehensive_dataset_final[, final_columns]

# Step 14: Save complete dataset
write.csv(comprehensive_dataset_final, "comprehensive_transcript_dataset.csv", row.names = FALSE)

# Step 15: Create filtered versions for paper
# Correlated + Discordant only
paper_dataset <- comprehensive_dataset_final %>%
  filter(is_correlated == TRUE | is_discordant == TRUE)

write.csv(paper_dataset, "final_paper_lists/paper_transcript_dataset.csv", row.names = FALSE)

# Step 16: Summary statistics
cat("=== Dataset Summary ===\n")
cat("Total transcripts:", nrow(comprehensive_dataset_final), "\n")
cat("Discordant transcripts:", sum(comprehensive_dataset_final$is_discordant, na.rm = TRUE), "\n") 
cat("Correlated transcripts:", sum(comprehensive_dataset_final$is_correlated, na.rm = TRUE), "\n")
cat("Protein coding:", sum(comprehensive_dataset_final$transcript_type == "protein_coding", na.rm = TRUE), "\n")
cat("Pass filtering:", sum(comprehensive_dataset_final$passes_filtering, na.rm = TRUE), "\n")
cat("Paper dataset size:", nrow(paper_dataset), "\n")

print("Comprehensive transcript dataset created successfully!")