library(dplyr)
library(tidyr)
library(ggplot2)

# Load CCLE expression data (using existing RDS files)
expression_gene <- readRDS("data/CCLE_gene_disease_Figure1B_new.rds")
transcript_expression <- readRDS("data/Figure_1C_CCLE_transcript_new.rds")

# Load Supplemental Table (comprehensive data set)
comprehensive_transcript_dataset <- read.csv("data/comprehensive_transcript_dataset.csv")

### CCLE GENE EXPRESSION DENSITY PLOTS ###

gene_list_from_transcripts <- comprehensive_transcript_dataset %>%
  pull(gene_id) %>%
  unique()

# Convert to data.table for speed
dt <- as.data.table(expression_gene)
gene_columns <- names(dt)[4:ncol(dt)]
clean_gene_names <- gsub(" \\(\\d+\\)$", "", gene_columns)

# Filter columns and melt
keep_cols <- gene_columns[clean_gene_names %in% gene_list_from_transcripts]
dt_subset <- dt[, c(names(dt)[1:3], keep_cols), with = FALSE]

long_expr <- melt(dt_subset, 
                 id.vars = names(dt)[1:3],
                 variable.name = "Gene", 
                 value.name = "Expression")
long_expr$Gene <- gsub(" \\(\\d+\\)$", "", long_expr$Gene)

long_expr_filtered <- long_expr %>%
  filter(Gene %in% gene_list_from_transcripts)

# Step 3: Calculate mean expression per DepmapModelType x Gene
summarized_expr <- long_expr_filtered %>%
  dplyr::group_by(DepmapModelType, Gene) %>%
  dplyr::summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Rank DepmapModelTypes per gene (lower rank = higher expression)
ranked_expr <- summarized_expr %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 5: Extract SKCM ranks only and add gene group classification
ranked_expr_joined <- ranked_expr %>%
  left_join(
    comprehensive_transcript_dataset %>%
      dplyr::select(gene_id, is_protein_coding_filtered, is_discordant, is_correlated) %>%
      distinct(gene_id, .keep_all = TRUE),
    by = c("Gene" = "gene_id")
  )

skcm_ranks <- ranked_expr_joined %>%
  filter(DepmapModelType == "SKCM") %>%
  mutate(GeneGroup = case_when(
    is_protein_coding_filtered ~ "Protein_Coding",
    is_discordant ~ "Discordant",
    is_correlated ~ "Correlated",
    TRUE ~ "Other"
  )) %>%
  filter(GeneGroup != "Other") %>%
  dplyr::select(Gene, MeanRank, GeneGroup)

# Step 6: Create density plots for each group - CCLE Gene
# Discordant genes
skcm_ranks_discordant <- skcm_ranks %>% filter(GeneGroup == "Discordant")
ggplot(skcm_ranks_discordant, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "blue", fill = "blue") +
  labs(
    title = "Density Plot of SKCM Ranks - Discordant Genes",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")

# Correlated genes  
skcm_ranks_correlated <- skcm_ranks %>% filter(GeneGroup == "Correlated")
ggplot(skcm_ranks_correlated, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "blue", fill = "blue") +
  labs(
    title = "Density Plot of SKCM Ranks - Correlated Genes",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")

# Protein coding genes
skcm_ranks_protein <- skcm_ranks %>% filter(GeneGroup == "Protein_Coding")
ggplot(skcm_ranks_protein, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "blue", fill = "blue") +
  labs(
    title = "Density Plot of SKCM Ranks - Protein Coding Genes",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")


### CCLE TRANSCRIPT EXPRESSION DENSITY PLOTS ###

# Step 1: Reshape transcript_expression to long format
long_expr_transcript <- transcript_expression %>%
  pivot_longer(cols = 4:61, names_to = "TranscriptID", values_to = "Expression")

# Step 2: Filter for your transcript lists
# Get transcript lists from comprehensive dataset
discordant_transcripts <- comprehensive_transcript_dataset %>%
  filter(is_discordant == TRUE) %>%
  pull(transcript_id)

correlated_transcripts <- comprehensive_transcript_dataset %>%
  filter(is_correlated == TRUE) %>%
  pull(transcript_id)

protein_coding_transcripts <- comprehensive_transcript_dataset %>%
  filter(is_protein_coding_filtered == TRUE) %>%
  pull(transcript_id)

# Then filter using these vectors
long_expr_transcript_filtered <- long_expr_transcript %>%
  filter(TranscriptID %in% c(discordant_transcripts,
                             correlated_transcripts,
                             protein_coding_transcripts))

# Step 3: Calculate mean expression per DepmapModelType x Transcript
summarized_expr_transcript <- long_expr_transcript_filtered %>%
  dplyr::group_by(DepmapModelType, TranscriptID) %>%
  dplyr::summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Rank DepmapModelTypes per transcript (lower rank = higher expression)
ranked_expr_transcript <- summarized_expr_transcript %>%
  dplyr::group_by(TranscriptID) %>%
  mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 5: Extract SKCM ranks only and add transcript group classification
skcm_ranks_transcript <- ranked_expr_transcript %>%
  filter(DepmapModelType == "SKCM") %>%
  mutate(TranscriptGroup = case_when(
    TranscriptID %in% transcript_unique$transcript_id ~ "Discordant",
    TranscriptID %in% transcript_correlation_all$transcript_id ~ "Correlated",
    TranscriptID %in% protein_coding$transcript_id ~ "Protein_Coding",
    TRUE ~ "Other"
  )) %>%
  filter(TranscriptGroup != "Other") %>%
  dplyr::select(TranscriptID, MeanRank, TranscriptGroup)

# Step 6: Create density plots for each group - CCLE Transcript
# Discordant transcripts
skcm_ranks_discordant_transcript <- skcm_ranks_transcript %>% filter(TranscriptGroup == "Discordant")
ggplot(skcm_ranks_discordant_transcript, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "red", fill = "red") +
  labs(
    title = "Density Plot of SKCM Ranks - Discordant Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")

# Correlated transcripts
skcm_ranks_correlated_transcript <- skcm_ranks_transcript %>% filter(TranscriptGroup == "Correlated")
ggplot(skcm_ranks_correlated_transcript, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "red", fill = "red") +
  labs(
    title = "Density Plot of SKCM Ranks - Correlated Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")

# Protein coding transcripts
skcm_ranks_protein_transcript <- skcm_ranks_transcript %>% filter(TranscriptGroup == "Protein_Coding")
ggplot(skcm_ranks_protein_transcript, aes(x = MeanRank)) +
  geom_density(alpha = 0.5, color = "red", fill = "red") +
  labs(
    title = "Density Plot of SKCM Ranks - Protein Coding Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +
  theme_minimal() +
  theme(legend.position = "none")
