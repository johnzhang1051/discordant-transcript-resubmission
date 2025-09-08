library(data.table)

# Import the large CSV using fread
Transcript_expression <- fread("Data/Transcript_expression_melanoma.csv")
correlated <- read.csv("Data/final_MITF_correlated.csv")

# Ensure both are data.tables
Transcript_expression <- as.data.table(Transcript_expression)
correlated <- as.data.table(correlated)
setnames(correlated, "transcript_ID", "transcript_id")

# Step 1: Remove gene name row
Transcript_expression_clean <- Transcript_expression[-1, ]

# Step 2: Keep only the Sample_ID column + transcripts in correlated
transcript_cols_to_keep <- intersect(names(Transcript_expression_clean), correlated$transcript_id)
Transcript_expression_subset <- Transcript_expression_clean[, c("Sample_ID", transcript_cols_to_keep), with = FALSE]

# Step 3: Convert transcript columns to numeric
Transcript_expression_subset[, (2:ncol(Transcript_expression_subset)) := 
                               lapply(.SD, as.numeric), .SDcols = 2:ncol(Transcript_expression_subset)]

# Step 4: Transpose (samples = rows → transcripts = rows)
expr_t <- transpose(Transcript_expression_subset[, -1])  # remove Sample_ID
colnames(expr_t) <- Transcript_expression_subset$Sample_ID  # set sample names
expr_t[, transcript_id := transcript_cols_to_keep]  # add transcript_id

# Step 5: TRUE/FALSE for >10 in ≥25% of samples
n_samples <- ncol(expr_t) - 1
expr_t[, pass_25pct := rowSums(.SD > 10, na.rm = TRUE) >= ceiling(0.25 * n_samples), .SDcols = 1:n_samples]

# Step 6: Merge with correlated
filtered_correlated <- merge(correlated, expr_t[, .(transcript_id, pass_25pct)], by = "transcript_id", all.x = TRUE)

# Step 7: Keep only passing
correlated_filtered <- filtered_correlated[pass_25pct == TRUE]



discordant_final_list <- read.csv("final_paper_LISTS/discordant_RESUBMISSION.csv")
library(dplyr)
correlated_PAPER <- correlated_filtered %>%
  filter(!(transcript_id %in% correlated_final_list$transcript_id))



  library(dplyr)

# Merge by transcript_id
transcript_type <- read.csv("Data/transcripttype.csv")

correlated_protein_coding <- correlated_PAPER %>%
  inner_join(transcript_type, by = "transcript_id")
correlated_protein_coding <- correlated_protein_coding %>%
  filter(transcript_type.x == "protein_coding")

if(nrow(correlated_PAPER) == nrow(correlated_protein_coding)) {
  write.csv(correlated_PAPER, "final_paper_LISTS/correlated_RESUBMISSION.csv", row.names = FALSE)
  print("Success")
} else {
  print("Not all transcripts are protein coding. Something weird has happened")
}

