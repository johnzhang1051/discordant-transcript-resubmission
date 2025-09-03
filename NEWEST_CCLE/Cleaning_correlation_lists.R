Pearson_Spearman_gene <- read.csv("PearsonSpearman_CCLE/ENST00000394351_gene_correlations.csv")
colnames(Pearson_Spearman_gene)[2:3] <- c("Spearman_gene", "Pearson_gene")


Spearman_transcript <- read.csv("PearsonSpearman_CCLE/cor_trans_spearman_protein_coding_clean.csv")
colnames(Spearman_transcript)[1:2] <- c("transcript_ID", "Spearman_transcript")

Pearson_transcript <- read.csv("PearsonSpearman_CCLE/cor_trans_pearson_protein_coding_clean.csv")
colnames(Pearson_transcript)[1:2] <- c("transcript_ID", "Pearson_transcript")

### get transcript to gene map
# Load biomaRt
library(biomaRt)
# Check available datasets to verify name
mart <- useEnsembl(biomart = "ensembl")
# List datasets if needed
# listDatasets(mart)
# Set human dataset
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# Download mapping of transcript ID to gene symbol
transcript_gene_map <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  mart = ensembl
)
# Optional: rename columns
colnames(transcript_gene_map) <- c("transcript_ID", "Gene")
# Save to CSV if needed
write.csv(transcript_gene_map, "PearsonSpearman_CCLE_clean/human_transcript_gene_map.csv", row.names = FALSE)





library(dplyr)

# Add Gene to Pearson_transcript
Pearson_transcript <- Pearson_transcript %>%
  left_join(transcript_gene_map, by = c("transcript_ID" = "transcript_ID"))

# Add Gene to Spearman_transcript
Spearman_transcript <- Spearman_transcript %>%
  left_join(transcript_gene_map, by = c("transcript_ID" = "transcript_ID"))

##
library(dplyr)
# Ensure dplyr is loaded
library(dplyr)

# Add Pearson_gene to Pearson_transcript by Gene
Pearson_transcript <- Pearson_transcript %>%
  dplyr::left_join(
    dplyr::select(Pearson_Spearman_gene, Gene, Pearson_gene),
    by = "Gene"
  )

# Add Spearman_transcript to Spearman_transcript by transcript_ID
Spearman_transcript <- Spearman_transcript %>%
  dplyr::left_join(
    dplyr::select(Pearson_Spearman_gene, Gene, Spearman_gene),
    by = "Gene"
  )


## remove NA values
library(dplyr)
library(tidyr)
# Remove rows with any NA values
Pearson_transcript_clean <- Pearson_transcript %>%
  drop_na()

Spearman_transcript_clean <- Spearman_transcript %>%
  drop_na()

# Save to CSV
write.csv(Pearson_transcript_clean, "PearsonSpearman_CCLE_clean/Pearson_transcript_clean.csv", row.names = FALSE)
write.csv(Spearman_transcript_clean, "PearsonSpearman_CCLE_clean/Spearman_transcript_clean.csv", row.names = FALSE)
