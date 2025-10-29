library(data.table)


## import and 
Melanoma_ID <- read.csv("Original_Data/Melanoma_Sample_ID.csv")
Transcript_expression <- fread("Original_Data/OmicsExpressionTranscriptsExpectedCountProfile.csv")
Gene_expression <- fread("Original_Data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
colnames(Gene_expression)[1] <- "Sample_ID"
colnames(Transcript_expression)[1] <- "Profile_ID"
Profile_to_sample <- read.csv("Original_Data/Profile_Sample_ID.csv")

## replace profile_ID for Sample_ID
library(dplyr)
# Assuming the structure:
# Transcript_expression: first column is Profile_ID, rest are expression values
# Profile_to_sample: columns Profile_ID and Sample_ID
# Step 1: Rename first column in Transcript_expression to "Profile_ID" (if not already)
colnames(Transcript_expression)[1] <- "Profile_ID"
# Step 2: Join with Profile_to_sample to get Sample_ID
Transcript_expression <- left_join(Transcript_expression, Profile_to_sample, by = "Profile_ID")
# Step 3: Move Sample_ID to the first column
Transcript_expression <- Transcript_expression %>%
  dplyr::select(Sample_ID, everything(), -Profile_ID)
# Step 4 (Optional): Rename Sample_ID to be the new first column name
colnames(Transcript_expression)[1] <- "Sample_ID"


# Ensure column is named consistently


# Subset for melanoma cell lines only
Gene_expression_melanoma <- Gene_expression[Gene_expression$Sample_ID %in% Melanoma_ID$Sample_ID, ]
Transcript_expression_melanoma <- Transcript_expression[Transcript_expression$Sample_ID %in% Melanoma_ID$Sample_ID, ]

# Subset gene expression to only samples in transcript data
Gene_expression_melanoma <- Gene_expression_melanoma %>%
  filter(Sample_ID %in% Transcript_expression_melanoma$Sample_ID)




## Clean up Transcript_expression_melanoma
library(stringr)


# Extract current column names
old_colnames <- colnames(Transcript_expression_melanoma)
first_col <- old_colnames[1]  # Typically "Sample_ID" or "Profile_ID"

# Parse gene and transcript IDs from remaining column names
split_result <- str_match(old_colnames[-1], "^(.+) \\((ENST[0-9]+)\\)$")

# Error check
if (nrow(split_result) != length(old_colnames[-1])) {
  stop("Regex match failed on some column names. Check format like 'GENE (ENST000...)'")
}

gene_names <- split_result[, 2]
tx_ids <- split_result[, 3]

# Rename columns to transcript IDs
colnames(Transcript_expression_melanoma) <- c(first_col, tx_ids)

# Convert gene names into a 1-row data frame with correct column names
gene_row_df <- as.data.frame(t(c(first_col, gene_names)), stringsAsFactors = FALSE)
colnames(gene_row_df) <- colnames(Transcript_expression_melanoma)

# Convert rest of data to character (if not already) and rbind
Transcript_expression_melanoma[] <- lapply(Transcript_expression_melanoma, as.character)
Transcript_expression_melanoma <- rbind(gene_row_df, Transcript_expression_melanoma)


### clean up gene expression
# Keep original first column name
first_col <- colnames(Gene_expression_melanoma)[1]

# Extract and clean remaining column names
gene_cols <- colnames(Gene_expression_melanoma)[-1]

# Remove " (number)" part using regex
clean_gene_names <- sub(" \\([0-9]+\\)$", "", gene_cols)

# Apply cleaned names to data frame
colnames(Gene_expression_melanoma) <- c(first_col, clean_gene_names)

## Save_Cleaned_Data_sets
# Step 1: Save Gene_expression_melanoma
write.csv(Gene_expression_melanoma, "Cleaned_Data/Gene_expression_melanoma.csv", row.names = FALSE)

# Step 2: Save Transcript_expression_melanoma
write.csv(Transcript_expression_melanoma, "Cleaned_Data/Transcript_expression_melanoma.csv", row.names = FALSE)

# Step 3: Create log2(x + 1) normalized version of Transcript_expression_melanoma

# Remove the first row (gene names) and convert values to numeric for log2
Transcript_expression_numeric <- Transcript_expression_melanoma[-1, ]
Transcript_expression_numeric[, -1] <- lapply(Transcript_expression_numeric[, -1], as.numeric)

# Apply log2(x + 1)
Transcript_expression_log2 <- Transcript_expression_numeric
Transcript_expression_log2[, -1] <- log2(Transcript_expression_numeric[, -1] + 1)

# Save normalized version
write.csv(Transcript_expression_log2, "Cleaned_Data/Transcript_expression_melanoma_log2.csv", row.names = FALSE)
