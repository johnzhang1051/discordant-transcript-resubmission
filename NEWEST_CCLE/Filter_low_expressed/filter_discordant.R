library(data.table)
library(conflicted)
conflicts_prefer(data.table::transpose)

# Import the large CSV using fread
Transcript_expression <- fread("Data/Transcript_expression_melanoma.csv")

discordant <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")

# Ensure both are data.tables
Transcript_expression <- as.data.table(Transcript_expression)
discordant <- as.data.table(discordant)

# Step 1: Remove gene name row
Transcript_expression_clean <- Transcript_expression[-1, ]

# Step 2: Keep only the Sample_ID column + transcripts in discordant
transcript_cols_to_keep <- intersect(names(Transcript_expression_clean), discordant$transcript_id)
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

# Step 6: Merge with discordant
filtered_discordant <- merge(discordant, expr_t[, .(transcript_id, pass_25pct)], by = "transcript_id", all.x = TRUE)

# Step 7: Keep only passing
discordant_filtered <- filtered_discordant[pass_25pct == TRUE]

write.csv(discordant_filtered, "final_paper_lists/discordant_RESUBMISSION.csv", row.names = FALSE)
