expression_gene <- readRDS("CCLE_gene_disease_Figure1B.rds")
transcript_expression <- readRDS("Figure_1C_CCLE_transcript.rds")


### CCLE MEDIAN EXPRESSION
library(tidyverse)

# Step 1: Select DepmapModelType and the expression columns
df_expr <- expression_gene[, c(3, 4:56)]

# Step 2: Calculate median expression per sample (row-wise) across gene columns only
df_expr <- df_expr %>%
  rowwise() %>%
  dplyr::mutate(median_expr = median(c_across(-DepmapModelType), na.rm = TRUE)) %>%
  ungroup()

# Step 3: Create a plotting data frame
plot_df <- df_expr %>%
  dplyr::select(DepmapModelType, median_expr) %>%
  dplyr::mutate(color_group = ifelse(DepmapModelType == "SKCM", "SKCM", "Other"))

# Step 4: Plot
ggplot(plot_df, aes(x = color_group, y = median_expr, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("SKCM" = "blue", "Other" = "gray70")) +
  labs(
    title = "Median Gene Expression in SKCM vs Other Tumor Types",
    x = "Tumor Type",
    y = "Median Gene Expression"
  ) +
  ylim(0, 8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


### CCLE transcript
library(tidyverse)

# Step 1: Select DepmapModelType and expression columns
df_expr <- transcript_expression[, c(3, 4:61)]

# Step 2: Calculate median expression per sample across genes
df_expr <- df_expr %>%
  rowwise() %>%
  dplyr::mutate(median_expr = median(c_across(-DepmapModelType), na.rm = TRUE)) %>%
  ungroup()

# Step 3: Prepare for plotting
plot_df <- df_expr %>%
  dplyr::select(DepmapModelType, median_expr) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "SKCM", "Other"))

# Step 4: Plot
ggplot(plot_df, aes(x = color_group, y = median_expr, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("SKCM" = "red", "Other" = "gray70")) +
  labs(
    title = "Median Transcript Expression in SKCM vs Other Tumor Types",
    x = "Tumor Type",
    y = "Median Transcript Expression"
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



### MEAN and MEDIAN ranked order CCLE gene

library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Reshape expression_gene to long format
long_expr <- expression_gene %>%
  pivot_longer(cols = 4:56, names_to = "TranscriptID", values_to = "Expression")

# Step 2: Calculate mean or median expression per DepmapModelType x Transcript
summarized_expr <- long_expr %>%
  group_by(DepmapModelType, TranscriptID) %>%
  dplyr::summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    MedianExpr = median(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Rank DepmapModelTypes per transcript (lower rank = higher expression)
ranked_expr <- summarized_expr %>%
  dplyr::group_by(TranscriptID) %>%
  dplyr::mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average"),
    MedianRank = rank(-MedianExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 4: Extract SKCM ranks only
skcm_ranks <- ranked_expr %>%
  dplyr::filter(DepmapModelType == "SKCM") %>%
  dplyr::select(TranscriptID, MeanRank, MedianRank)

# Step 5: Plot rank distribution
ggplot(skcm_ranks, aes(x = reorder(TranscriptID, MeanRank), y = MeanRank)) +
  geom_point(color = "blue", alpha = 0.6) +
  labs(
    title = "SKCM Rank Order Across Transcripts (Mean Expression)",
    x = "Transcript (ordered by SKCM rank)",
    y = "SKCM Rank (across DepmapModelTypes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Step 6: Summary statistic
mean_skcm_rank <- mean(skcm_ranks$MeanRank)
median_skcm_rank <- median(skcm_ranks$MedianRank)

cat("Average SKCM rank (mean-based):", mean_skcm_rank, "\n")
cat("Median SKCM rank (median-based):", median_skcm_rank, "\n")

### nice box plot
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming `skcm_ranks` is already created as before:

# Convert to long format for ggplot
skcm_ranks_long <- skcm_ranks %>%
  pivot_longer(cols = c(MeanRank, MedianRank), names_to = "RankType", values_to = "Rank")
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(fill = "blue", alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  theme_minimal()

ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_density(alpha = 0.5, color = "blue") +
  scale_fill_manual(values = c("blue", "blue")) +
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis scale without cutting off data
  labs(
    title = "Density Plot of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")




ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("blue", "blue")) +
  labs(
    title = "Histogram of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


### box plot with blue points
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(
    fill = "blue", color = "blue", width = 0.6, alpha = 0.4, outlier.shape = NA
  ) +
  geom_point(
    position = position_nudge(x = 0),  # no horizontal jitter
    color = "blue", alpha = 0.5, size = 1
  ) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  coord_cartesian(ylim = c(1, 16)) +
  theme_minimal()


### MEAN and MEDIAN ranked order transcript   CCLE

library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Reshape expression_gene to long format
long_expr <- transcript_expression %>%
  pivot_longer(cols = 4:61, names_to = "TranscriptID", values_to = "Expression")

# Step 2: Calculate mean or median expression per DepmapModelType x Transcript
summarized_expr <- long_expr %>%
  group_by(DepmapModelType, TranscriptID) %>%
  summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    MedianExpr = median(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Rank DepmapModelTypes per transcript (lower rank = higher expression)
ranked_expr <- summarized_expr %>%
  group_by(TranscriptID) %>%
  mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average"),
    MedianRank = rank(-MedianExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 4: Extract SKCM ranks only
skcm_ranks <- ranked_expr %>%
  filter(DepmapModelType == "SKCM") %>%
  select(TranscriptID, MeanRank, MedianRank)

# Step 5: Plot rank distribution
ggplot(skcm_ranks, aes(x = reorder(TranscriptID, MeanRank), y = MeanRank)) +
  geom_point(color = "red", alpha = 0.6) +
  labs(
    title = "SKCM Rank Order Across Transcripts (Mean Expression)",
    x = "Transcript (ordered by SKCM rank)",
    y = "SKCM Rank (across DepmapModelTypes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Step 6: Summary statistic
mean_skcm_rank <- mean(skcm_ranks$MeanRank)
median_skcm_rank <- median(skcm_ranks$MedianRank)

cat("Average SKCM rank (mean-based):", mean_skcm_rank, "\n")
cat("Median SKCM rank (median-based):", median_skcm_rank, "\n")

### nice box plot
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming `skcm_ranks` is already created as before:

# Convert to long format for ggplot
skcm_ranks_long <- skcm_ranks %>%
  pivot_longer(cols = c(MeanRank, MedianRank), names_to = "RankType", values_to = "Rank")
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(fill = "blue", alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  theme_minimal()

ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_density(alpha = 0.5, color = "red") +
  scale_fill_manual(values = c("red", "red")) +
  labs(
    title = "Density Plot of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +  # More stable than xlim()
  theme_minimal() +
  theme(legend.position = "none")


ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("red", "red")) +
  labs(
    title = "Histogram of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

### box plot with blue points
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(
    fill = "red", color = "red", width = 0.6, alpha = 0.4, outlier.shape = NA
  ) +
  geom_point(
    position = position_nudge(x = 0),  # no horizontal jitter
    color = "red", alpha = 0.5, size = 1
  ) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  coord_cartesian(ylim = c(1, 16)) +
  theme_minimal()



### import TCGA transcript expression
TCGA_transcript_expression <- readRDS("TCGA_expr_clean.rds")
# 1. Extract the primary disease row (row 1) from expression columns
primary_disease <- TCGA_transcript_expression[1, 3:ncol(TCGA_transcript_expression)]
# 2. Remove the first row and the "Gene" column
expr_data <- TCGA_transcript_expression[-1, -2]  # drop first row and column 2 ("Gene")
# 3. Transpose the expression matrix
expr_transposed <- as.data.frame(t(expr_data[, -1]))  # remove transcript_id before transpose
# 4. Set column names as transcript IDs
colnames(expr_transposed) <- expr_data[[1]]  # first column has transcript IDs
# 5. Add sample ID as first column
expr_transposed <- cbind(Sample_ID = colnames(expr_data)[-1], expr_transposed)
# 6. Add primary disease as second column
expr_transposed <- cbind(Primary_Disease = as.character(primary_disease), expr_transposed)
# 7. Reorder columns: Sample_ID first, Primary_Disease second, then transcript expression
expr_transposed <- expr_transposed[, c("Sample_ID", "Primary_Disease", setdiff(names(expr_transposed), c("Sample_ID", "Primary_Disease")))]

### MEAN and MEDIAN ranked order transcript   TCGA

library(dplyr)
library(tidyr)
library(ggplot2)
transcript_expression <- expr_transposed
# Step 1: Reshape expression_gene to long format
long_expr <- transcript_expression %>%
  pivot_longer(cols = 3:55, names_to = "TranscriptID", values_to = "Expression")

# Step 2: Calculate mean or median expression per Primary Disease x Transcript
summarized_expr <- long_expr %>%
  mutate(Expression = as.numeric(Expression)) %>%
  group_by(Primary_Disease, TranscriptID) %>%
  summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    MedianExpr = median(Expression, na.rm = TRUE),
    .groups = "drop"
  )


# Step 3: Rank DepmapModelTypes per transcript (lower rank = higher expression)
ranked_expr <- summarized_expr %>%
  group_by(TranscriptID) %>%
  mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average"),
    MedianRank = rank(-MedianExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 4: Extract SKCM ranks only
skcm_ranks <- ranked_expr %>%
  filter(Primary_Disease == "SKCM") %>%
  select(TranscriptID, MeanRank, MedianRank)

# Step 5: Plot rank distribution
ggplot(skcm_ranks, aes(x = reorder(TranscriptID, MeanRank), y = MeanRank)) +
  geom_point(color = "red", alpha = 0.6) +
  labs(
    title = "SKCM Rank Order Across Transcripts (Mean Expression)",
    x = "Transcript (ordered by SKCM rank)",
    y = "SKCM Rank (across DepmapModelTypes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Step 6: Summary statistic
mean_skcm_rank <- mean(skcm_ranks$MeanRank)
median_skcm_rank <- median(skcm_ranks$MedianRank)

cat("Average SKCM rank (mean-based):", mean_skcm_rank, "\n")
cat("Median SKCM rank (median-based):", median_skcm_rank, "\n")

### nice box plot
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming `skcm_ranks` is already created as before:

# Convert to long format for ggplot
skcm_ranks_long <- skcm_ranks %>%
  pivot_longer(cols = c(MeanRank, MedianRank), names_to = "RankType", values_to = "Rank")
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(fill = "red", alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  theme_minimal()

ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_density(alpha = 0.5, color = "red") +
  scale_fill_manual(values = c("red", "red")) +
  labs(
    title = "Density Plot of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(1, 16)) +  # More stable than xlim()
  theme_minimal() +
  theme(legend.position = "none")


ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("red", "red")) +
  labs(
    title = "Histogram of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

### box plot with blue points
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(
    fill = "red", color = "red", width = 0.6, alpha = 0.4, outlier.shape = NA
  ) +
  geom_point(
    position = position_nudge(x = 0),  # no horizontal jitter
    color = "red", alpha = 0.5, size = 1
  ) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  coord_cartesian(ylim = c(1, 16)) +
  theme_minimal()




## TCGA Gene Expression

### import TCGA transcript expression
TCGA_transcript_expression <- readRDS("TCGA_gene_expr_clean.rds")
# 1. Extract the primary disease row (row 1) from expression columns
primary_disease <- TCGA_transcript_expression[1, 3:ncol(TCGA_transcript_expression)]
# 2. Remove the first row and the "Gene" column
expr_data <- TCGA_transcript_expression[-1, -2]  # drop first row and column 2 ("Gene")
# 3. Transpose the expression matrix
expr_transposed <- as.data.frame(t(expr_data[, -1]))  # remove transcript_id before transpose
# 4. Set column names as transcript IDs
colnames(expr_transposed) <- expr_data[[1]]  # first column has transcript IDs
# 5. Add sample ID as first column
expr_transposed <- cbind(Sample_ID = colnames(expr_data)[-1], expr_transposed)
# 6. Add primary disease as second column
expr_transposed <- cbind(Primary_Disease = as.character(primary_disease), expr_transposed)
# 7. Reorder columns: Sample_ID first, Primary_Disease second, then transcript expression
expr_transposed <- expr_transposed[, c("Sample_ID", "Primary_Disease", setdiff(names(expr_transposed), c("Sample_ID", "Primary_Disease")))]

### MEAN and MEDIAN ranked order Gene   TCGA

library(dplyr)
library(tidyr)
library(ggplot2)
transcript_expression <- expr_transposed
# Step 1: Reshape expression_gene to long format
long_expr <- transcript_expression %>%
  pivot_longer(cols = 3:55, names_to = "TranscriptID", values_to = "Expression")

# Step 2: Calculate mean or median expression per Primary Disease x Transcript
summarized_expr <- long_expr %>%
  mutate(Expression = as.numeric(Expression)) %>%
  group_by(Primary_Disease, TranscriptID) %>%
  summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    MedianExpr = median(Expression, na.rm = TRUE),
    .groups = "drop"
  )


# Step 3: Rank DepmapModelTypes per transcript (lower rank = higher expression)
ranked_expr <- summarized_expr %>%
  group_by(TranscriptID) %>%
  mutate(
    MeanRank = rank(-MeanExpr, ties.method = "average"),
    MedianRank = rank(-MedianExpr, ties.method = "average")
  ) %>%
  ungroup()

# Step 4: Extract SKCM ranks only
skcm_ranks <- ranked_expr %>%
  filter(Primary_Disease == "SKCM") %>%
  select(TranscriptID, MeanRank, MedianRank)

# Step 5: Plot rank distribution
ggplot(skcm_ranks, aes(x = reorder(TranscriptID, MeanRank), y = MeanRank)) +
  geom_point(color = "red", alpha = 0.6) +
  labs(
    title = "SKCM Rank Order Across Transcripts (Mean Expression)",
    x = "Transcript (ordered by SKCM rank)",
    y = "SKCM Rank (across DepmapModelTypes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Step 6: Summary statistic
mean_skcm_rank <- mean(skcm_ranks$MeanRank)
median_skcm_rank <- median(skcm_ranks$MedianRank)

cat("Average SKCM rank (mean-based):", mean_skcm_rank, "\n")
cat("Median SKCM rank (median-based):", median_skcm_rank, "\n")

### nice box plot
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming `skcm_ranks` is already created as before:

# Convert to long format for ggplot
skcm_ranks_long <- skcm_ranks %>%
  pivot_longer(cols = c(MeanRank, MedianRank), names_to = "RankType", values_to = "Rank")
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(fill = "blue", alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  theme_minimal()

ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_density(alpha = 0.5, color = "blue") +
  scale_fill_manual(values = c("blue", "blue")) +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 0.25)) +  # Set both axes
  labs(
    title = "Density Plot of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



ggplot(skcm_ranks_long, aes(x = Rank, fill = RankType)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("blue", "blue")) +
  labs(
    title = "Histogram of SKCM Ranks Across Transcripts",
    x = "SKCM Rank (Lower = More Highly Expressed)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

### box plot with blue points
ggplot(skcm_ranks_long, aes(x = RankType, y = Rank)) +
  geom_boxplot(
    fill = "blue", color = "blue", width = 0.6, alpha = 0.4, outlier.shape = NA
  ) +
  geom_point(
    position = position_nudge(x = 0),  # no horizontal jitter
    color = "blue", alpha = 0.5, size = 1
  ) +
  labs(
    title = "Distribution of SKCM Ranks Across Transcripts",
    x = "Rank Type",
    y = "SKCM Rank (Lower = More Highly Expressed)"
  ) +
  coord_cartesian(ylim = c(1, 16)) +
  theme_minimal()



