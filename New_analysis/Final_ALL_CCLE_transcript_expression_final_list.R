# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(cowplot)

# === 1. LOAD DATA ===
expr_df <- readRDS("Data/CCLE_transcript_expression.rds")
# Extract column names (skip first column, 'X')
col_names <- colnames(expr_df)[-1]

# Split on double period
split_names <- strsplit(col_names, "\\.\\.")

# Extract gene and transcript
gene_names <- sapply(split_names, `[`, 1)
transcript_ids <- sapply(split_names, `[`, 2)

# Create header rows
header1 <- c(NA, gene_names)
header2 <- c(NA, transcript_ids)

# Grab the data body (including sample names in first column)
data_body <- expr_df

# Combine all as a data frame
full_df <- rbind(header1, header2, as.matrix(data_body))
full_df <- as.data.frame(full_df, stringsAsFactors = FALSE)

# Optional: set meaningful column names (e.g., V1, V2, ...)
colnames(full_df) <- c("Sample", paste0("V", 2:ncol(full_df)))

# Save as RDS
saveRDS(full_df, file = "CCLE_transcript_needs_cleaned.rds")

### START HERE
expr_temp <- readRDS("Data/CCLE_transcript_needs_cleaned.rds")
colnames(expr_temp)[1] <- "ModelID"
model_ID <-  read.csv("Data/cell_line_disease.csv")

# Step 1: Extract gene names and transcript IDs
gene_names <- as.character(expr_temp[1, ])
transcript_ids <- as.character(expr_temp[2, ])

# Step 2: Set transcript IDs as column names
colnames(expr_temp) <- transcript_ids
colnames(expr_temp)[1] <- "ModelID"  # label first column properly

# Step 3: Remove the first two rows (we've saved them)
expr_temp_data <- expr_temp[-c(1, 2), ]

# Step 4: Convert to data frame and fix ModelID column
expr_temp_data <- as.data.frame(expr_temp_data)
expr_temp_data$ModelID <- as.character(expr_temp_data$ModelID)
colnames(expr_temp_data) <- transcript_ids
colnames(expr_temp_data)[1] <- "ModelID"

# Step 5: Convert expression columns to numeric
expr_temp_data[ , -1] <- lapply(expr_temp_data[ , -1], as.numeric)

# Step 6: Merge with model_id
merged_expr <- merge(model_ID, expr_temp_data, by = "ModelID")
colnames(merged_expr)[4:198840] <- sub("\\.$", "", colnames(merged_expr)[4:198840])


# Step 7: Create transcript → gene map
transcript_gene_map <- data.frame(
  TranscriptID = transcript_ids[-1],  # exclude "ModelID"
  Gene = gene_names[-1],
  stringsAsFactors = FALSE
)






TCGA_modeltype <- read.csv("Data/TCGA_DepmapModelType.txt")
TCGA_modeltype_unique <- unique(trimws(as.character(TCGA_modeltype$DepmapModelType)))
# Clean up the column for matching
merged_expr$DepmapModelType <- trimws(as.character(merged_expr$DepmapModelType))
# Filter to keep only matching model types
filtered_CCLE_transcript <- merged_expr[merged_expr$DepmapModelType %in% TCGA_modeltype_unique, ]
final_list <- read.csv("Data/final_transcript_list.txt")


# Step 1: Extract transcript IDs from final_list and clean
final_transcripts <- as.character(final_list$transcript_id)
final_transcripts <- sub("\\.$", "", trimws(final_transcripts))  # remove trailing dots, spaces

# Step 2: Clean column names in filtered_CCLE_transcript
colnames(filtered_CCLE_transcript) <- sub("\\.$", "", trimws(colnames(filtered_CCLE_transcript)))

# Step 3: Find matching transcript columns
matching_cols <- which(colnames(filtered_CCLE_transcript) %in% final_transcripts)

# Step 4: Keep metadata (cols 1–3) and matching expression columns
filtered_df <- filtered_CCLE_transcript[, c(1:3, matching_cols[matching_cols > 3])]


library(ggplot2)
library(dplyr)
library(forcats)

# Assign color: SKCM = pink, others = gray
rank_summary_df$Color <- ifelse(rank_summary_df$DepmapModelType == "SKCM", "pink", "gray")

# --- Summarize mean & median across all transcripts ---
summary_means <- rank_summary_df %>%
  group_by(DepmapModelType, Color) %>%
  summarise(Value = mean(mean_expr, na.rm = TRUE)) %>%
  mutate(Metric = "Mean")

summary_medians <- rank_summary_df %>%
  group_by(DepmapModelType, Color) %>%
  summarise(Value = mean(median_expr, na.rm = TRUE)) %>%
  mutate(Metric = "Median")

# --- Combine into one data frame ---
summary_combined <- bind_rows(summary_means, summary_medians)

# --- Reorder by Value (use mean for ordering) ---
ordering <- summary_means %>%
  arrange(Value) %>%
  pull(DepmapModelType)

summary_combined$DepmapModelType <- factor(summary_combined$DepmapModelType, levels = ordering)

# --- Plot both metrics side-by-side ---
combined_plot <- ggplot(summary_combined, aes(x = DepmapModelType, y = Value, fill = Color)) +
  geom_col(position = position_dodge(width = 0.8), aes(group = Metric), width = 0.7) +
  facet_wrap(~Metric) +
  scale_fill_identity() +
  theme_bw() +
  labs(title = "Ranked Mean & Median Transcript Expression by Model Type",
       x = "DepMap Model Type", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# --- Save to PDF ---
ggsave("RankSummaries/Combined_Transcript_Expression_By_ModelType.pdf", combined_plot, width = 10, height = 6)


## rank based summarization
library(ggplot2)
library(dplyr)
library(forcats)

# Step 1: Compute rank per transcript
ranked_df <- rank_summary_df %>%
  group_by(Gene) %>%
  mutate(Rank = rank(-mean_expr, ties.method = "average"))  # "-" for descending order (1 = highest expression)

# Step 2: Average rank per DepmapModelType
average_ranks <- ranked_df %>%
  group_by(DepmapModelType) %>%
  summarise(AverageRank = mean(Rank)) %>%
  mutate(Color = ifelse(DepmapModelType == "SKCM", "pink", "gray"))

# Step 3: Order by average rank (low to high = higher expression)
average_ranks$DepmapModelType <- fct_reorder(average_ranks$DepmapModelType, average_ranks$AverageRank)

# Step 4: Plot
rank_plot <- ggplot(average_ranks, aes(x = DepmapModelType, y = AverageRank, fill = Color)) +
  geom_col() +
  scale_fill_identity() +
  theme_bw() +
  labs(title = "Average Rank of Transcript Expression by Model Type",
       x = "DepMap Model Type", y = "Average Rank (Lower = Higher Expression)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 5: Save plot and data
ggsave("RankSummaries/Average_Transcript_Rank_By_ModelType.pdf", rank_plot, width = 9, height = 6)
write.csv(average_ranks, "RankSummaries/Average_Transcript_Rank_By_ModelType.csv", row.names = FALSE)

### ABR PLOT FOR PAPER
### prepare filtered DF as RDS
filtered_df <- filtered_CCLE_transcript[, c(1:3, matching_cols[matching_cols > 3])]
saveRDS(filtered_df, "Figure_1C_CCLE_transcript.rds")

### import df and plot ABR
filtered_df_for_ABR_CCLE_transcript <- readRDS("Figure_1C_CCLE_transcript.rds")


# Prepare the data
abr_df <- filtered_df_for_ABR_CCLE_transcript %>%
  select(DepmapModelType = 3, ENST00000544583) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000544583, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ENST00000544583, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
