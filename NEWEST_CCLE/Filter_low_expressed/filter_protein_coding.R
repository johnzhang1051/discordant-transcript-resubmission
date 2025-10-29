library(data.table)

# Import the large CSV using fread
Transcript_expression <- fread("Data/Transcript_expression_melanoma.csv")

transcript_type <- read.csv("Data/transcripttype.csv")

protein_coding <- transcript_type %>% filter(transcript_type == "protein_coding")


# Ensure both are data.tables
Transcript_expression <- as.data.table(Transcript_expression)
protein_coding <- as.data.table(protein_coding)

# Step 1: Remove gene name row
Transcript_expression_clean <- Transcript_expression[-1, ]

# Step 2: Keep only the Sample_ID column + transcripts in protein_coding
transcript_cols_to_keep <- intersect(names(Transcript_expression_clean), protein_coding$transcript_id)
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

# Step 6: Merge with protein_coding
filtered_protein_coding <- merge(protein_coding, expr_t[, .(transcript_id, pass_25pct)], by = "transcript_id", all.x = TRUE)

# Step 7: Keep only passing
protein_coding_filtered <- filtered_protein_coding[pass_25pct == TRUE]

write.csv(protein_coding_filtered, "final_paper_lists/protein_coding_RESUBMISSION.csv", row.names = FALSE)
