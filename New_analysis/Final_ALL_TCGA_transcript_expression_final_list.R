# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(cowplot)

# === 1. LOAD DATA ===
expr_df <- readRDS("Data/TCGA_transcript_cleaned.rds")
transcript_ids <- read.csv("Data/final_list_transcript_id.txt")
TCGA_phenotype <- read.csv("Data/TCGA_phenotype.csv")
Gene_transcript <- read.csv("Data/final_transcript_list.csv")

 

# === 2. FILTER & JOIN ===
colnames(expr_df)[1] <- "transcript_id"
expr_subset <- inner_join(transcript_ids, expr_df, by = "transcript_id")
expr_annotated <- left_join(expr_subset, Gene_transcript, by = "transcript_id")
expr_annotated <- expr_annotated %>% select(transcript_id, Gene, everything())

# === 3. ADD PRIMARY DISEASE AS FIRST ROW ===
sample_cols <- setdiff(colnames(expr_annotated), c("transcript_id", "Gene"))
phenotype_matched <- TCGA_phenotype[TCGA_phenotype$sample %in% sample_cols, ]
primary_disease_row <- phenotype_matched$primary_disease[match(sample_cols, phenotype_matched$sample)]
new_row <- c("primary_disease", NA, primary_disease_row)
new_row_df <- as.data.frame(t(new_row), stringsAsFactors = FALSE)
colnames(new_row_df) <- colnames(expr_annotated)
expr_annotated_with_disease <- rbind(new_row_df, expr_annotated)

# === 4. CONVERT TO LONG FORMAT ===
expr_data <- expr_annotated_with_disease[-1, ]
expr_data <- expr_data %>% mutate(across(-c(transcript_id, Gene), as.numeric))
long_expr <- expr_data %>%
  pivot_longer(cols = -c(transcript_id, Gene), names_to = "sample", values_to = "expression") %>%
  left_join(TCGA_phenotype, by = c("sample" = "sample")) %>%
  filter(!is.na(expression), !is.na(primary_disease))

# === 5. DEFINE COLOR PALETTE ===
unique_diseases <- unique(long_expr$primary_disease)
disease_palette <- setNames(rep("gray70", length(unique_diseases)), unique_diseases)
disease_palette["SKCM"] <- "#ff1493"  # Bright pink
disease_palette["UVM"]  <- "#8a2be2"  # Bright purple
other_diseases <- setdiff(unique_diseases, c("SKCM", "UVM"))
safe_colors <- brewer.pal(min(length(other_diseases), 8), "Set2")
disease_palette[other_diseases] <- safe_colors

# === 6. MULTI-PAGE PDF: EXPRESSION PLOTS ===
pdf("transcript_expression_by_disease.pdf", width = 10, height = 6)
transcript_ids <- unique(long_expr$transcript_id)

for (tid in transcript_ids) {
  df <- long_expr %>% filter(transcript_id == tid)
  if (nrow(df) == 0) next
  
  df <- df %>%
    group_by(primary_disease) %>%
    mutate(median_expr = median(expression, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(primary_disease = fct_reorder(primary_disease, median_expr))
  
  gene_name <- unique(df$Gene)
  
  p <- ggplot(df, aes(x = primary_disease, y = expression, fill = primary_disease)) +
    geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5) +
    scale_fill_manual(values = disease_palette) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = "Primary Disease",
      y = "Expression (TPM)",
      fill = "Disease Type",
      title = paste0("Gene: ", gene_name, "\nTranscript ID: ", tid)
    )
  print(p)
}
dev.off()

# === 7. CALCULATE MEAN + MEDIAN RANK PER DISEASE ===
ranked_df <- long_expr %>%
  group_by(transcript_id, primary_disease) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(transcript_id) %>%
  mutate(rank = rank(-mean_expr, ties.method = "average")) %>%
  ungroup()

rank_summary <- ranked_df %>%
  group_by(primary_disease) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value") %>%
  mutate(highlight = ifelse(primary_disease %in% c("SKCM", "UVM"), "highlight", "other"),
         fill_key = paste(statistic, highlight, sep = "."))

color_map <- c(
  "mean_rank.highlight" = "#ff1493",
  "median_rank.highlight" = "#8a2be2",
  "mean_rank.other" = "gray70",
  "median_rank.other" = "gray50"
)

# === 8. PLOT RANK SUMMARY AND SAVE ===
p_rank <- ggplot(rank_long, aes(x = reorder(primary_disease, rank_value), y = rank_value, fill = fill_key)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Mean and Median Rank of Expression by Disease",
    x = "Primary Disease",
    y = "Rank (lower = higher expression)",
    fill = "Statistic"
  ) +
  scale_fill_manual(
    values = color_map,
    labels = c(
      "mean_rank.highlight" = "Mean Rank (SKCM/UVM)",
      "mean_rank.other" = "Mean Rank (Other)",
      "median_rank.highlight" = "Median Rank (SKCM/UVM)",
      "median_rank.other" = "Median Rank (Other)"
    )
  )

ggsave("expression_rank_summary.pdf", plot = p_rank, width = 10, height = 6)

# Save rank summary table to CSV
write.csv(rank_summary, "expression_rank_summary_table.csv", row.names = FALSE)

# Optional: also save as RDS for reloading in R
saveRDS(rank_summary, "expression_rank_summary_table.rds")


### save to RDS for script for paper
saveRDS(expr_annotated_with_disease, "TCGA_transcript_disease_Figure1B.rds")


### ABR expression
TCGA_transcript_expr_annotated_with_disease <- readRDS("TCGA_transcript_disease_Figure1B.rds")
# Vector of primary diseases to keep
keep_diseases <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA",
                   "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD",
                   "UCS", "UCEC", "HNSC", "SKCM")

# Extract disease labels from row 1
sample_disease <- as.character(TCGA_transcript_expr_annotated_with_disease[1, 3:ncol(TCGA_transcript_expr_annotated_with_disease)])
names(sample_disease) <- colnames(TCGA_transcript_expr_annotated_with_disease)[3:ncol(TCGA_transcript_expr_annotated_with_disease)]

# Subset to samples with desired primary diseases
samples_to_keep <- names(sample_disease)[sample_disease %in% keep_diseases]

# Clean expression dataframe
expr_clean <- TCGA_transcript_expr_annotated_with_disease[1:54, ]

# Filter for ABR
abr_expr <- expr_clean %>%
  filter(Gene == "ABR") %>%
  select(transcript_id, Gene, all_of(samples_to_keep))

# Pivot to long format
abr_long <- abr_expr %>%
  pivot_longer(
    cols = -c(transcript_id, Gene),
    names_to = "Sample",
    values_to = "Expression"
  )

# Add disease label and color group
abr_long <- abr_long %>%
  mutate(Disease = sample_disease[Sample],
         color_group = ifelse(Disease == "SKCM", "Melanoma", "Other"),
         Expression = as.numeric(Expression))

# Reorder diseases by median ABR expression (highest on right)
abr_long <- abr_long %>%
  mutate(Disease = fct_reorder(Disease, Expression, .fun = median))

# Plot
ggplot(abr_long, aes(x = Disease, y = Expression, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "Tumor Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

