✅ Assumptions:
  You have a target variable y: a numeric vector of length n_samples.

You have:
  
  Transcript: a data frame where rows = samples, column 1 is Sample_ID, remaining columns are transcript expression with column names as transcript_ID.

Gene: same format, but for gene expression.

You have a transcript ↔ gene mapping: Pearson_transcript_gene with transcript_ID, Gene.

library(data.table)
Pearson_transcript_gene <- read.csv("CoCor_Data_sets/human_transcript_gene_map.csv")
Transcript <- fread("CoCor_Data_sets/Transcript_expression_melanoma_log2.csv")
Gene <- fread("CoCor_Data_sets/Gene_expression_melanoma.csv")
Transcript <- Transcript[Transcript$Sample_ID != "ACH-000931", ]
Gene <- Gene[Gene$Sample_ID != "ACH-000931", ]

# Load required packages
library(dplyr)
library(purrr)
library(tibble)
library(cocor)

# Ensure correct rownames on expression matrices
trans_expr <- Transcript %>% column_to_rownames("Sample_ID")
gene_expr  <- Gene %>% column_to_rownames("Sample_ID")
common_samples <- intersect(rownames(trans_expr), rownames(gene_expr))
trans_expr <- trans_expr[common_samples, ]
gene_expr <- gene_expr[common_samples, ]
all(rownames(trans_expr) == rownames(gene_expr))  # must be TRUE


# Define y as outcome (expression of the target transcript)
y <- trans_expr[["ENST00000394351"]]

# Filter valid pairs: transcript in trans_expr AND gene in gene_expr
valid_pairs <- Pearson_transcript_gene %>%
  filter(transcript_ID %in% colnames(trans_expr),
         Gene %in% colnames(gene_expr))

# Function to run cocor for one transcript-gene pair
run_cocor_test <- function(trans_id, gene_id, y_vec) {
  x1 <- trans_expr[[trans_id]]
  x2 <- gene_expr[[gene_id]]
  
  complete_idx <- complete.cases(x1, x2, y_vec)
  x1 <- x1[complete_idx]
  x2 <- x2[complete_idx]
  y  <- y_vec[complete_idx]
  
  if (length(y) < 10) {
    return(data.frame(
      transcript_ID = trans_id,
      Gene = gene_id,
      p_value = NA,
      n = length(y),
      r_transcript = NA,
      r_gene = NA,
      r_transcript_gene = NA
    ))
  }
  
  r.jk <- cor(x1, y, method = "pearson")
  r.jh <- cor(x2, y, method = "pearson")
  r.kh <- cor(x1, x2, method = "pearson")
  
  
  test <- tryCatch({
    cocor.dep.groups.overlap(r.jk = r.jk, r.jh = r.jh, r.kh = r.kh, n = length(y))
  }, error = function(e) NULL)
  
  p_val <- tryCatch({
    if (!is.null(test)) {
      as.numeric(test@steiger1980$p.value)
    } else {
      NA
    }
  }, error = function(e) NA)
  
  
  
  return(data.frame(
    transcript_ID = trans_id,
    Gene = gene_id,
    p_value = p_val,
    n = length(y),
    r_transcript = r.jk,
    r_gene = r.jh,
    r_transcript_gene = r.kh
  ))
}

# Apply function to all valid transcript-gene pairs
results <- purrr::pmap_dfr(
  valid_pairs %>% select(transcript_ID, Gene),
  function(transcript_ID, Gene) run_cocor_test(transcript_ID, Gene, y)
)

# Optional: Filter significant results
results_filtered <- results %>% filter(!is.na(p_value) & p_value < 0.05)

# Save results
write.csv(results, "cocor_all_results.csv", row.names = FALSE)
write.csv(results_filtered, "cocor_significant_results.csv", row.names = FALSE)

### add FDR correction
results <- results %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

results_fdr <- results %>% filter(!is.na(FDR) & FDR < 0.05)

# Save results
write.csv(results_fdr, "cocor_fdr_significant_results.csv", row.names = FALSE)

# Subset where both correlations are ≥ 0.5
results_fdr_r0.5_both <- results_fdr %>%
  filter(!is.na(r_transcript), !is.na(r_gene),
         r_transcript >= 0.5, r_gene <= 0.5)

# Save to CSV
write.csv(results_fdr_r0.5_both, "cocor_fdr_rtranscript_rgene_0.5.csv", row.names = FALSE)

### same code for Spearman
run_cocor_test_spearman <- function(trans_id, gene_id, y_vec) {
  x1 <- trans_expr[[trans_id]]
  x2 <- gene_expr[[gene_id]]
  
  complete_idx <- complete.cases(x1, x2, y_vec)
  x1 <- x1[complete_idx]
  x2 <- x2[complete_idx]
  y  <- y_vec[complete_idx]
  
  if (length(y) < 10) {
    return(data.frame(
      transcript_ID = trans_id,
      Gene = gene_id,
      p_value = NA,
      n = length(y),
      r_transcript = NA,
      r_gene = NA,
      r_transcript_gene = NA
    ))
  }
  
  r.jk <- cor(x1, y, method = "spearman")
  r.jh <- cor(x2, y, method = "spearman")
  r.kh <- cor(x1, x2, method = "spearman")
  
  test <- tryCatch({
    cocor.dep.groups.overlap(r.jk = r.jk, r.jh = r.jh, r.kh = r.kh, n = length(y))
  }, error = function(e) NULL)
  
  p_val <- tryCatch({
    if (!is.null(test)) {
      as.numeric(test@steiger1980$p.value)
    } else {
      NA
    }
  }, error = function(e) NA)
  
  return(data.frame(
    transcript_ID = trans_id,
    Gene = gene_id,
    p_value = p_val,
    n = length(y),
    r_transcript = r.jk,
    r_gene = r.jh,
    r_transcript_gene = r.kh
  ))
}

results_spearman <- purrr::pmap_dfr(
  valid_pairs %>% select(transcript_ID, Gene),
  function(transcript_ID, Gene) run_cocor_test_spearman(transcript_ID, Gene, y)
)

# Add FDR column
results_spearman <- results_spearman %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

# Filter FDR-significant results
results_spearman_fdr <- results_spearman %>%
  filter(!is.na(FDR) & FDR < 0.05)

# Subset 1: r_transcript ≥ 0.5
spearman_r_transcript_0.5 <- results_spearman_fdr %>%
  filter(!is.na(r_transcript), r_transcript >= 0.5)

# Subset 2: r_transcript ≥ 0.5 and r_gene ≥ 0.5
spearman_r_transcript_gene_0.5 <- results_spearman_fdr %>%
  filter(!is.na(r_transcript), !is.na(r_gene),
         r_transcript >= 0.5, r_gene >= 0.5)

# Subset 3: r_transcript ≥ 0.5 and r_gene < 0.5
spearman_r_transcript_0.5_gene_lt_0.5 <- results_spearman_fdr %>%
  filter(!is.na(r_transcript), !is.na(r_gene),
         r_transcript >= 0.5, r_gene < 0.5)

write.csv(results_spearman, "cocor_spearman_all_results.csv", row.names = FALSE)
write.csv(results_spearman_fdr, "cocor_spearman_fdr_results.csv", row.names = FALSE)
write.csv(spearman_r_transcript_0.5, "cocor_spearman_rtranscript_0.5.csv", row.names = FALSE)
write.csv(spearman_r_transcript_0.5_gene_lt_0.5, "cocor_spearman_rtranscript_0.5_rgene_lt_0.5.csv", row.names = FALSE)

merged_transcripts <- bind_rows(
  spearman_r_transcript_0.5_gene_lt_0.5,
  results_fdr_r0.5_both
) %>%
  distinct(transcript_ID, .keep_all = TRUE)

write.csv(merged_transcripts, "cocor_CCLE_final_list.csv")


transcript_type <- read.csv("CoCor_Data_sets/transcripttype.csv")

# Ensure dplyr is loaded
library(dplyr)

# Merge and filter for protein_coding transcripts
CoCor_final_protein_coding <- merged_transcripts %>%
  left_join(transcript_type, by = "transcript_ID") %>%
  filter(transcript_type == "protein_coding")

