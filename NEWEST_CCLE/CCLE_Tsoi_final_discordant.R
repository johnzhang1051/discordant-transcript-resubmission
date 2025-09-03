library(tximport)
library(readr)
library(dplyr)

# Step 1: Build transcript-to-gene mapping (if not already done)
# NOTE: If you already have this as a .tsv or prebuilt data.frame, you can skip this
gtf_file <- "gencode.v44.annotation.gtf"
tx2gene <- rtracklayer::import(gtf_file) %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  dplyr::select(transcript_id = transcript_id, gene_id = gene_id) %>%
  distinct()

# Step 1: List only valid subfolders (exclude hidden ones)
sample_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\./\\.", sample_dirs)]  # remove .hidden folders

# Step 2: Keep only those that contain 'abundance.tsv'
files <- file.path(sample_dirs, "abundance.tsv")
files <- files[file.exists(files)]  # remove paths where the file does not exist
names(files) <- basename(dirname(files))

# Check final list
print(files)
library(rhdf5)
# Step 3: Run tximport
txi_gene <- tximport(files,
                     type = "kallisto",
                     tx2gene = tx2gene,
                     countsFromAbundance = "no",
                     ignoreAfterBar = TRUE)


# Step 5: Run tximport for transcript-level counts
txi_tx <- tximport(files,
                   type = "kallisto",
                   txOut = TRUE,
                   countsFromAbundance = "no",
                   ignoreAfterBar = TRUE)

library(tibble)

# Step 6: Extract and save gene-level counts
gene_counts <- as.data.frame(txi_gene$counts) %>%
  rownames_to_column("gene_id")
write.csv(gene_counts, "Tsoi/kallisto_gene_counts.csv", row.names = FALSE)

# Step 7: Extract and save transcript-level counts
transcript_counts <- as.data.frame(txi_tx$counts) %>%
  rownames_to_column("transcript_id")
write.csv(transcript_counts, "Tsoi/kallisto_transcript_counts.csv", row.names = FALSE)

# Optional: Save TPMs too
gene_tpm <- as.data.frame(txi_gene$abundance) %>%
  rownames_to_column("gene_id")
gene_tpm$gene_id <- sub("\\..*", "", gene_tpm$gene_id)
write.csv(gene_tpm, "Tsoi/kallisto_gene_TPM.csv", row.names = FALSE)

transcript_tpm <- as.data.frame(txi_tx$abundance) %>%
  rownames_to_column("transcript_id")
# Trim version numbers (remove everything after the first ".")
transcript_tpm$transcript_id <- sub("\\..*", "", transcript_tpm$transcript_id)
write.csv(transcript_tpm, "Tsoi/kallisto_transcript_TPM.csv", row.names = FALSE)

####### START HERE JZ ######
setwd("/Users/johnz/Documents/GitFiles/discordant-transcript-resubmission/NEWEST_CCLE")
transcript_tpm <- read.csv("Tsoi/kallisto_transcript_TPM.csv")
transcript_tpm$transcript_id <- sub("\\..*", "", transcript_tpm$transcript_id)


### spearman and pearson
# Step 1: Prepare matrix (remove ID column and log2-transform)
expr_mat <- transcript_tpm
rownames(expr_mat) <- expr_mat$transcript_id
expr_mat <- expr_mat[, -1]  # remove transcript_id column
expr_mat_log <- log2(expr_mat + 1)

# Step 2: Get expression vector for reference transcript
ref_id <- "ENST00000394351"

if (!ref_id %in% rownames(expr_mat_log)) {
  stop("Reference transcript not found.")
}

ref_expr <- expr_mat_log[ref_id, ]
# Transpose so rows = samples, columns = transcripts
expr_mat_log_t <- t(expr_mat_log)

# Get the reference transcript expression vector
ref_expr <- expr_mat_log_t[, "ENST00000394351"]

# Compute correlations
pearson_r <- apply(expr_mat_log_t, 2, function(x) cor(x, ref_expr, method = "pearson"))
spearman_r <- apply(expr_mat_log_t, 2, function(x) cor(x, ref_expr, method = "spearman"))


# Step 4: Combine results into a data frame
cor_results <- data.frame(
  transcript_id = rownames(expr_mat_log),
  Pearson = pearson_r,
  Spearman = spearman_r
)

# Step 5: Save to CSV
write.csv(cor_results, "Tsoi/correlation_to_ENST00000394351.csv", row.names = FALSE)

# Step 6: Filter for Pearson ≥ 0.5 OR Spearman ≥ 0.5
cor_selected <- cor_results %>%
  dplyr::filter(Pearson >= 0.5 | Spearman >= 0.5)

# Step 7: Save selected transcripts
write.csv(cor_selected, "Tsoi/correlated_transcripts_r0.5.csv", row.names = FALSE)

# Optional: Extract just the list of transcript IDs
writeLines(cor_selected$transcript_id, "transcript_ids_r0.5.txt")

CCLE_discordant <- read.csv("CoCor_analysis/CoCor_final_list_protein_coding.csv")
colnames(CCLE_discordant)[2] <- "transcript_id"

library(dplyr)

# Make sure CCLE_discordant has a column named "transcript_id"
# If not, rename it:
# colnames(CCLE_discordant)[<column_number>] <- "transcript_id"

# Step 1: Find overlapping transcripts
overlap_transcripts <- inner_join(cor_selected, CCLE_discordant, by = "transcript_id")

# Step 2: Save the overlap
write.csv(overlap_transcripts, "Tsoi/CCLE_Tsoi_discordant_protein_coding.csv", row.names = FALSE)
