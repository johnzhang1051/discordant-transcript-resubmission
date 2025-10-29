library(dplyr)
library(conflicted)


conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

# Function to convert columns to numeric where possible
convert_to_numeric <- function(df) {
  for (col in names(df)) {
    # Try to convert to numeric
    numeric_col <- suppressWarnings(as.numeric(df[[col]]))
    
    # If conversion resulted in mostly non-NA values, keep it
    # (Allow up to 10% NAs from conversion)
    original_nas <- sum(is.na(df[[col]]))
    new_nas <- sum(is.na(numeric_col))
    
    if ((new_nas - original_nas) / nrow(df) < 0.1) {
      df[[col]] <- numeric_col
      cat("Converted column '", col, "' to numeric\n", sep = "")
    }
  }
  return(df)
}

# Comparison function with automatic conversion
compare_csv_summary <- function(file1, file2) {
  
  # Read the files
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  
  cat("=" , rep("=", 70), "\n", sep = "")
  cat("CSV FILES COMPARISON\n")
  cat("=" , rep("=", 70), "\n\n", sep = "")
  
  # Convert to numeric where possible
  cat("Converting File 1 columns to numeric...\n")
  df1 <- convert_to_numeric(df1)
  cat("\nConverting File 2 columns to numeric...\n")
  df2 <- convert_to_numeric(df2)
  cat("\n")
  
  # Basic dimensions
  cat("DIMENSIONS:\n")
  cat("  File 1:", nrow(df1), "rows x", ncol(df1), "columns\n")
  cat("  File 2:", nrow(df2), "rows x", ncol(df2), "columns\n\n")
  
  # Column names
  cat("COLUMNS:\n")
  cat("  File 1:", paste(colnames(df1), collapse = ", "), "\n")
  cat("  File 2:", paste(colnames(df2), collapse = ", "), "\n\n")
  
  # Only compare if same columns exist
  common_cols <- intersect(colnames(df1), colnames(df2))
  
  if (length(common_cols) > 0) {
    cat("SUMMARY STATISTICS FOR COMMON COLUMNS:\n\n")
    
    for (col in common_cols) {
      cat("Column:", col, "\n")
      
      if (is.numeric(df1[[col]]) && is.numeric(df2[[col]])) {
        cat("  File 1 - Mean:", round(mean(df1[[col]], na.rm = TRUE), 3),
            "| Median:", round(median(df1[[col]], na.rm = TRUE), 3),
            "| SD:", round(sd(df1[[col]], na.rm = TRUE), 3),
            "| Min:", round(min(df1[[col]], na.rm = TRUE), 3),
            "| Max:", round(max(df1[[col]], na.rm = TRUE), 3),
            "| NAs:", sum(is.na(df1[[col]])), "\n")
        cat("  File 2 - Mean:", round(mean(df2[[col]], na.rm = TRUE), 3),
            "| Median:", round(median(df2[[col]], na.rm = TRUE), 3),
            "| SD:", round(sd(df2[[col]], na.rm = TRUE), 3),
            "| Min:", round(min(df2[[col]], na.rm = TRUE), 3),
            "| Max:", round(max(df2[[col]], na.rm = TRUE), 3),
            "| NAs:", sum(is.na(df2[[col]])), "\n")
        
        # Show difference in means
        diff <- mean(df1[[col]], na.rm = TRUE) - mean(df2[[col]], na.rm = TRUE)
        if (abs(diff) > 0.001) {
          cat("  *** Mean difference:", round(diff, 3), "\n")
        }
      } else {
        cat("  File 1 - Unique values:", length(unique(df1[[col]])),
            "| NAs:", sum(is.na(df1[[col]])), "\n")
        cat("  File 2 - Unique values:", length(unique(df2[[col]])),
            "| NAs:", sum(is.na(df2[[col]])), "\n")
      }
      cat("\n")
    }
  }
  
  invisible(list(df1 = df1, df2 = df2))
}

# Usage
# result <- compare_csv_summary("file1.csv", "file2.csv")

# Read and convert
df1 <- read.csv("CoCor_analysis/CoCor_final_list_protein_coding_old.csv")
df2 <- read.csv("CoCor_analysis/CoCor_final_list_protein_coding.csv")



df1 <- convert_to_numeric(df1)
df2 <- convert_to_numeric(df2)

# Now compare
summary(df1)
summary(df2)


# Create the filtered discordant lists
df1_discordant <- df1 %>% 
  filter(r_transcript >= 0.5, r_gene < 0.5)

df2_discordant <- df2 %>% 
  filter(r_transcript >= 0.5, r_gene < 0.5)

# Identify lost transcripts
lost_transcript_ids <- setdiff(df1_discordant$transcript_ID, df2_discordant$transcript_ID)

cat("=== LOST TRANSCRIPTS ANALYSIS ===\n\n")
cat("Total lost:", length(lost_transcript_ids), "transcripts\n\n")

# Get the old values (from df1) for lost transcripts
lost_old_values <- df1 %>% 
  filter(transcript_ID %in% lost_transcript_ids)

cat("OLD VALUES (df1) - When they PASSED the filter:\n")
cat("--------------------------------------------------\n")
print(summary(lost_old_values[, c("r_transcript", "r_gene", "r_transcript_gene", "p_value")]))

# Get the new values (from df2) for lost transcripts, if they exist
lost_new_values <- df2 %>% 
  filter(transcript_ID %in% lost_transcript_ids)

if (nrow(lost_new_values) > 0) {
  cat("\n\nNEW VALUES (df2) - Why they FAILED the filter:\n")
  cat("--------------------------------------------------\n")
  print(summary(lost_new_values[, c("r_transcript", "r_gene", "r_transcript_gene", "p_value")]))
  
  # Check which criteria failed
  cat("\n\nFAILURE BREAKDOWN:\n")
  cat("Transcripts with r_transcript < 0.5:", sum(lost_new_values$r_transcript < 0.5, na.rm = TRUE), "\n")
  cat("Transcripts with r_gene >= 0.5:", sum(lost_new_values$r_gene >= 0.5, na.rm = TRUE), "\n")
  cat("Transcripts failing BOTH:", sum(lost_new_values$r_transcript < 0.5 & lost_new_values$r_gene >= 0.5, na.rm = TRUE), "\n")
  
  # Show the delta (change) in values
  comparison <- lost_old_values %>%
    select(transcript_ID, Gene, r_transcript_old = r_transcript, r_gene_old = r_gene) %>%
    inner_join(
      lost_new_values %>% select(transcript_ID, r_transcript_new = r_transcript, r_gene_new = r_gene),
      by = "transcript_ID"
    ) %>%
    mutate(
      delta_r_transcript = r_transcript_new - r_transcript_old,
      delta_r_gene = r_gene_new - r_gene_old
    )
  
  cat("\n\nCHANGE IN VALUES (New - Old):\n")
  cat("--------------------------------------------------\n")
  cat("Delta r_transcript:\n")
  print(summary(comparison$delta_r_transcript))
  cat("\nDelta r_gene:\n")
  print(summary(comparison$delta_r_gene))
  
  # Show worst offenders
  cat("\n\nTOP 10 TRANSCRIPTS WITH BIGGEST r_gene INCREASE:\n")
  cat("--------------------------------------------------\n")
  top_r_gene_increase <- comparison %>%
    arrange(desc(delta_r_gene)) %>%
    head(10) %>%
    select(transcript_ID, Gene, r_gene_old, r_gene_new, delta_r_gene)
  print(top_r_gene_increase)
  
  cat("\n\nTOP 10 TRANSCRIPTS WITH BIGGEST r_transcript DECREASE:\n")
  cat("--------------------------------------------------\n")
  top_r_transcript_decrease <- comparison %>%
    arrange(delta_r_transcript) %>%
    head(10) %>%
    select(transcript_ID, Gene, r_transcript_old, r_transcript_new, delta_r_transcript)
  print(top_r_transcript_decrease)
  
} else {
  cat("\n\nALL LOST TRANSCRIPTS WERE COMPLETELY REMOVED FROM df2!\n")
  cat("These transcripts don't exist in the new dataset at all.\n")
  cat("This means they failed the >10 TPM filter in 25% of samples.\n")
}

# Compare retained vs lost transcripts in OLD data
retained_transcript_ids <- intersect(df1_discordant$transcript_ID, df2_discordant$transcript_ID)
retained_old_values <- df1 %>% 
  filter(transcript_ID %in% retained_transcript_ids)

cat("\n\n=== COMPARISON: RETAINED vs LOST (in OLD data) ===\n\n")

cat("RETAINED transcripts (n =", length(retained_transcript_ids), ") - OLD values:\n")
print(summary(retained_old_values[, c("r_transcript", "r_gene", "r_transcript_gene")]))

cat("\n\nLOST transcripts (n =", length(lost_transcript_ids), ") - OLD values:\n")
print(summary(lost_old_values[, c("r_transcript", "r_gene", "r_transcript_gene")]))


df_all_results <- read.csv("CoCor_analysis/cocor_all_results.csv")
# Get lost transcript IDs
df1_discordant <- df1 %>% 
  filter(r_transcript >= 0.5, r_gene < 0.5)

df2_discordant <- df2 %>% 
  filter(r_transcript >= 0.5, r_gene < 0.5)

lost_transcript_ids <- setdiff(df1_discordant$transcript_ID, df2_discordant$transcript_ID)

# Get old values for lost transcripts
lost_old <- df1 %>% 
  filter(transcript_ID %in% lost_transcript_ids) %>%
  select(transcript_ID, Gene, 
         r_transcript_old = r_transcript, 
         r_gene_old = r_gene,
         r_transcript_gene_old = r_transcript_gene,
         p_value_old = p_value)

# Get new values from df_all_results
lost_comparison <- lost_old %>%
  left_join(
    df_all_results %>% 
      select(transcript_ID = transcript_ID, 
             r_transcript_new = r_transcript, 
             r_gene_new = r_gene,
             r_transcript_gene_new = r_transcript_gene,
             p_value_new = p_value),
    by = "transcript_ID"
  ) %>%
  mutate(
    # Calculate changes
    delta_r_transcript = r_transcript_new - r_transcript_old,
    delta_r_gene = r_gene_new - r_gene_old,
    delta_r_transcript_gene = r_transcript_gene_new - r_transcript_gene_old,
    delta_p_value = p_value_new - p_value_old,
    
    # Percent change
    pct_change_r_transcript = 100 * (r_transcript_new - r_transcript_old) / abs(r_transcript_old),
    pct_change_r_gene = 100 * (r_gene_new - r_gene_old) / abs(r_gene_old),
    
    # Label borderline status
    r_transcript_borderline = ifelse(r_transcript_old >= 0.5 & r_transcript_old < 0.55, "YES", "NO"),
    r_gene_borderline = ifelse(r_gene_old >= 0.4 & r_gene_old < 0.5, "YES", "NO"),
    borderline_category = case_when(
      r_transcript_old < 0.55 & r_gene_old >= 0.4 ~ "BOTH",
      r_transcript_old < 0.55 ~ "r_transcript only",
      r_gene_old >= 0.4 ~ "r_gene only",
      TRUE ~ "Neither (robust)"
    ),
    
    # Determine failure reason
    failed_r_transcript = r_transcript_new < 0.5,
    failed_r_gene = r_gene_new >= 0.5,
    failure_reason = case_when(
      is.na(r_transcript_new) ~ "Missing from new data",
      failed_r_transcript & failed_r_gene ~ "Both failed",
      failed_r_transcript ~ "r_transcript < 0.5",
      failed_r_gene ~ "r_gene >= 0.5",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(transcript_ID, Gene, borderline_category, failure_reason,
         r_transcript_old, r_transcript_new, delta_r_transcript, pct_change_r_transcript,
         r_gene_old, r_gene_new, delta_r_gene, pct_change_r_gene,
         r_transcript_gene_old, r_transcript_gene_new, delta_r_transcript_gene,
         p_value_old, p_value_new, delta_p_value,
         everything())

# View the comparison
View(lost_comparison)
