library(dplyr)
library(data.table)
library(tidyr)

# Load transcript lists
discordant_list <- read.csv("data/final_paper_lists/discordant_RESUBMISSION.csv")
correlated_list <- read.csv("data/final_paper_lists/correlated_RESUBMISSION.csv")


# Load correlation data
ccle_pearson <- read.csv("data/CoCor_analysis/cocor_all_results.csv")
ccle_spearman <- read.csv("data/CoCor_analysis/cocor_spearman_all_results.csv")
tsoi_transcript_correlations <- read.csv("data/Tsoi/correlation_to_ENST00000394351.csv")

# Tsoi gene correlations?

# Start with all unique transcripts
all_transcripts <- unique(c(
  ccle_pearson$transcript_ID,
  ccle_spearman$transcript_ID,
  tsoi_transcript_correlations$transcript_id
))


supplemental_table <- data.frame(
  transcript_id = all_transcripts,
  stringsAsFactors = FALSE
)

# Add CCLE correlations
ccle_corr_data <- ccle_pearson %>%
  select(transcript_id = transcript_ID,
         ccle_transcript_pearson = r_transcript,
         ccle_gene_pearson = r_gene) %>%
  left_join(
    ccle_spearman %>%
      select(transcript_id = transcript_ID,
             ccle_transcript_spearman = r_transcript,
             ccle_gene_spearman = r_gene),
    by = "transcript_id"
  )

supplemental_table <- left_join(supplemental_table, ccle_corr_data, by = "transcript_id")

# Add Tsoi transcript-level correlations
tsoi_corr_data <- tsoi_transcript_correlations %>%
  select(transcript_id,
         tsoi_transcript_pearson = Pearson,
         tsoi_transcript_spearman = Spearman)

supplemental_table <- left_join(supplemental_table, tsoi_corr_data, by = "transcript_id")

# Add MITF Overexpression data
GSE163646_OE <- read.csv("data/mitf_oe/GSE_163646_OE_kallisto_transcript_TPM.csv")
PRJNA704810_OE <- read.csv("data/mitf_oe/PRJNA704810_OE_kallisto_transcript_TPM.csv")

merged_OE <- merge(GSE163646_OE, PRJNA704810_OE, by = "transcript_id")

# Define sample groups
prjna_con <- c("SRR13782518", "SRR13782519")
prjna_oe <- c("SRR13782520", "SRR13782521", "SRR13782522")
gse_oe <- c("SRR13282354", "SRR13282355", "SRR13282356")
gse_con <- c("SRR13282357", "SRR13282358", "SRR13282359")

# Calculate OE ratios
oe_ratios <- merged_OE %>%
  mutate(
    PRJNA_CON_mean = rowMeans(select(., all_of(prjna_con)), na.rm = TRUE),
    PRJNA_OE_mean = rowMeans(select(., all_of(prjna_oe)), na.rm = TRUE),
    GSE_CON_mean = rowMeans(select(., all_of(gse_con)), na.rm = TRUE),
    GSE_OE_mean = rowMeans(select(., all_of(gse_oe)), na.rm = TRUE),
    PRJNA704810_OE_MITF_CON = PRJNA_OE_mean / PRJNA_CON_mean,
    GSE163646_OE_MITF_CON = GSE_OE_mean / GSE_CON_mean
  ) %>%
  select(transcript_id, PRJNA704810_OE_MITF_CON, GSE163646_OE_MITF_CON) %>%
  filter(is.finite(PRJNA704810_OE_MITF_CON), is.finite(GSE163646_OE_MITF_CON))

supplemental_table <- left_join(supplemental_table, oe_ratios, by = "transcript_id")

# Add MITF Knockdown data
PRJEB30337 <- read.csv("data/mitf_oe/PRJEB30337_TPM.csv")
GSE163646_KD <- read.csv("data/mitf_oe/GSE_163646_kallisto_transcript_TPM.csv")
Henja <- read.csv("data/mitf_oe/Henja_TPM.csv")
GSE283655 <- read.csv("data/mitf_oe/GSE283655_kallisto_transcript_TPM.csv")

# Clean transcript IDs
Henja$transcript_id <- sub("\\.\\d+$", "", Henja$transcript_id)
PRJEB30337$transcript_id <- sub("\\.\\d+$", "", PRJEB30337$transcript_id)

# Merge knockdown datasets
merged_siMITF <- Reduce(function(x, y) merge(x, y, by = "transcript_id", all = TRUE),
                        list(PRJEB30337, GSE163646_KD, Henja, GSE283655))

# Define knockdown sample groups
PRJEB30337_con <- c("ERR3013912", "ERR3013913")
PRJEB30337_siMITF <- c("ERR3013908", "ERR3013909")
GSE163646_CON <- c("SRR13282348", "SRR13282350")
GSE163646_siMITF <- c("SRR13282352", "SRR13282353")
GSE283655_con <- c("SRR31631254", "SRR31631258", "SRR31631262")
GSE283655_siMITF <- c("SRR31631252", "SRR31631256", "SRR31631260")
Henja_con <- c("SRR7346993", "SRR7346994")
Henja_siMITF <- c("SRR7346991", "SRR7346992")

# Calculate knockdown ratios (siMITF/siCON)
kd_ratios <- merged_siMITF %>%
  mutate(
    PRJEB30337_CON_mean = rowMeans(select(., all_of(PRJEB30337_con)), na.rm = TRUE),
    PRJEB30337_siMITF_mean = rowMeans(select(., all_of(PRJEB30337_siMITF)), na.rm = TRUE),
    GSE163646_CON_mean = rowMeans(select(., all_of(GSE163646_CON)), na.rm = TRUE),
    GSE163646_siMITF_mean = rowMeans(select(., all_of(GSE163646_siMITF)), na.rm = TRUE),
    GSE283655_CON_mean = rowMeans(select(., all_of(GSE283655_con)), na.rm = TRUE),
    GSE283655_siMITF_mean = rowMeans(select(., all_of(GSE283655_siMITF)), na.rm = TRUE),
    Henja_CON_mean = rowMeans(select(., all_of(Henja_con)), na.rm = TRUE),
    Henja_siMITF_mean = rowMeans(select(., all_of(Henja_siMITF)), na.rm = TRUE),
    PRJEB30337_siCON_siMITF = PRJEB30337_siMITF_mean / PRJEB30337_CON_mean,
    GSE163646_siCON_siMITF = GSE163646_siMITF_mean / GSE163646_CON_mean,
    GSE283655_siCON_siMITF = GSE283655_siMITF_mean / GSE283655_CON_mean,
    GSE115845_siCON_siMITF = Henja_siMITF_mean / Henja_CON_mean
  ) %>%
  select(transcript_id, PRJEB30337_siCON_siMITF, GSE163646_siCON_siMITF, 
         GSE283655_siCON_siMITF, GSE115845_siCON_siMITF)

supplemental_table <- left_join(supplemental_table, kd_ratios, by = "transcript_id")

# Add ChIP-seq peak data
Kenny <- read.csv("data/chip/Kenny.csv")
Laurette <- read.csv("data/chip/Laurette.csv")
Louph <- read.csv("data/chip/Louphrasitthiphol.csv")

# Standardize transcript IDs
standardize_transcripts <- function(df) {
  df %>%
    rename(transcript_id = transcriptId) %>%
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))
}

Kenny <- standardize_transcripts(Kenny)
Laurette <- standardize_transcripts(Laurette)
Louph <- standardize_transcripts(Louph)

# Filter for promoter peaks
filter_promoters <- function(df) {
  df %>%
    filter(annotation == "Promoter (<=1kb)") %>%
    select(transcript_id, annotation) %>%
    distinct()
}

Kenny_promoters <- filter_promoters(Kenny)
Laurette_promoters <- filter_promoters(Laurette)
Louph_promoters <- filter_promoters(Louph)

# Add peak information
supplemental_table <- supplemental_table %>%
  mutate(
    has_peak_GSE153020_Kenny = transcript_id %in% Kenny_promoters$transcript_id,
    has_peak_GSE61965_Laurette = transcript_id %in% Laurette_promoters$transcript_id,
    has_peak_GSE77437_Louph = transcript_id %in% Louph_promoters$transcript_id
  )

# Add classification flags
supplemental_table <- supplemental_table %>%
  filter(is_correlated | is_discordant)

supplemental_table <- supplemental_table %>%
  mutate(
    is_correlated = transcript_id %in% correlated_list$transcript_id,
    is_discordant = transcript_id %in% discordant_list$transcript_id
  )

# Add gene ID
# Extract gene_id from correlated list
correlated_gene_mapping <- correlated_list %>%
  select(transcript_id, gene_id = Gene.x) %>%
  distinct()

# Extract gene_id from discordant list
discordant_gene_mapping <- discordant_list %>%
  select(transcript_id, gene_id = Gene) %>%
  distinct()

# Combine both mappings (discordant takes priority if there's overlap)
gene_mapping_combined <- bind_rows(
  correlated_gene_mapping,
  discordant_gene_mapping
) %>%
  distinct(transcript_id, .keep_all = TRUE)  # Keep first occurrence (or use group_by logic)

# Join gene_id into supplemental_table
# This will update existing gene_id or add it if missing
supplemental_table <- supplemental_table %>%
  left_join(gene_mapping_combined, by = "transcript_id")

# Load CCLE expression data
ccle_transcript_expression <- fread("data/Transcript_expression_melanoma_log2.csv")
ccle_gene_expression <- fread("data/Gene_expression_melanoma.csv")

# Load Tsoi expression data
Tsoi_trans_expr <- read.csv("data/Tsoi/kallisto_transcript_TPM.csv")
Tsoi_gene_expr <- read.csv("data/Tsoi/kallisto_gene_TPM.csv")

# Calculate CCLE median transcript expression for each transcript
ccle_median_trans <- ccle_transcript_expression %>%
  pivot_longer(-Sample_ID, names_to = "transcript_id", values_to = "expression") %>%
  group_by(transcript_id) %>%
  summarise(ccle_median_transcript_expression = median(expression, na.rm = TRUE), .groups = "drop")

# Calculate Tsoi median transcript expression for each transcript
# Log2 transform to match CCLE
Tsoi_trans_log <- Tsoi_trans_expr %>%
  pivot_longer(-transcript_id, names_to = "sample", values_to = "expression") %>%
  mutate(expression = log2(expression + 1)) %>%
  group_by(transcript_id) %>%
  summarise(tsoi_median_transcript_expression = median(expression, na.rm = TRUE), .groups = "drop")

# Calculate CCLE median gene expression for each gene
ccle_median_gene <- ccle_gene_expression %>%
  pivot_longer(-Sample_ID, names_to = "gene_id", values_to = "expression") %>%
  group_by(gene_id) %>%
  summarise(ccle_median_gene_expression = median(expression, na.rm = TRUE), .groups = "drop")

# Calculate Tsoi median gene expression for each gene
# Log2 transform to match CCLE
Tsoi_gene_log <- Tsoi_gene_expr %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  mutate(expression = log2(expression + 1)) %>%
  group_by(gene_id) %>%
  summarise(tsoi_median_gene_expression = median(expression, na.rm = TRUE), .groups = "drop")

# Join all median expression data to supplemental_table
supplemental_table <- supplemental_table %>%
  left_join(ccle_median_trans, by = "transcript_id") %>%
  left_join(Tsoi_trans_log, by = "transcript_id") %>%
  left_join(ccle_median_gene, by = "gene_id") %>%
  left_join(Tsoi_gene_log, by = "gene_id")

# Reorder columns to match requested order
column_order <- c(
  "transcript_id",
  "gene_id",
  "ccle_transcript_spearman",
  "ccle_transcript_pearson",
  "tsoi_transcript_spearman",
  "tsoi_transcript_pearson",
  "ccle_gene_spearman",
  "ccle_gene_pearson",
  "tsoi_gene_spearman",
  "tsoi_gene_pearson",
  "PRJNA704810_OE_MITF_CON",
  "GSE163646_OE_MITF_CON",
  "PRJEB30337_siCON_siMITF",
  "GSE163646_siCON_siMITF",
  "GSE283655_siCON_siMITF",
  "GSE115845_siCON_siMITF",
  "has_peak_GSE153020_Kenny",
  "has_peak_GSE61965_Laurette",
  "has_peak_GSE77437_Louph",
  "transcript_type",
  "is_correlated",
  "is_discordant"
)

# Select only columns that exist
existing_columns <- intersect(column_order, colnames(supplemental_table))
supplemental_table_final <- supplemental_table[, existing_columns]

# Save the table
write.csv(supplemental_table_final, "supplemental_table.csv", row.names = FALSE)

######################## Expression data ########################

# Load CCLE expression data
ccle_transcript_expression <- fread("data/Transcript_expression_melanoma_log2.csv")
ccle_gene_expression <- fread("data/Gene_expression_melanoma.csv")

# Load Tsoi expression data
Tsoi_trans_expr <- read.csv("data/Tsoi/kallisto_transcript_TPM.csv")
Tsoi_gene_expr <- read.csv("data/Tsoi/kallisto_gene_TPM.csv")

# Export expression by transcript, just get median MITF expression per transcript
tsoi_transcript_mitf_expression <- Tsoi_trans_expr %>%
  mutate(
    MITF_median_expression = median(as.numeric(mitf_raw[1, -1]), na.rm = TRUE),
    MITF_mean_expression = mean(as.numeric(mitf_raw[1, -1]), na.rm = TRUE)
  ) %>%
  select(transcript_id, MITF_median_expression, MITF_mean_expression)
