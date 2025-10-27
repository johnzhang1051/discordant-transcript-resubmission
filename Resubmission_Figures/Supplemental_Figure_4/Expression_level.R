library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load data
CCLE_trans_expr <- read_tsv("data/CCLE_Melanoma_transcript_expression.txt")
CCLE_trans_expr <- CCLE_trans_expr %>%
  mutate(transcript_id = sub("\\..*$", "", transcript_id))

CCLE_trans_expr_clean <- CCLE_trans_expr[,-1]

# Step 1: Calculate the row means (average counts per transcript)
CCLE_trans_expr_clean$mean_count <- rowMeans(CCLE_trans_expr_clean[, 2:64], na.rm = TRUE)

# Step 2: Log2 transform (mean + 1)
CCLE_trans_expr_clean$log2_mean_plus1 <- log2(CCLE_trans_expr_clean$mean_count + 1)

CCLE_trans_expr_log2 <- CCLE_trans_expr_clean[,-(2:65)]

# Load transcript groups
correlation <- read.csv("data/correlated_RESUBMISSION.csv")
final <- read.csv("data/discordant_RESUBMISSION.csv")
protein_coding <- read.csv("data/protein_coding_RESUBMISSION.csv")

# Get protein-coding IDs
protein_coding_ids <- unique(protein_coding$transcript_id)
correlation_protein_ids <- unique(correlation$transcript_id)
final_protein_ids <- unique(final$transcript_id)


# Assign group labels for all three groups (protein-coding only)
expr_with_group <- CCLE_trans_expr_log2 %>%
  filter(transcript_id %in% protein_coding_ids) %>%
  mutate(group = case_when(
    transcript_id %in% final_protein_ids ~ "Unique",
    transcript_id %in% correlation_protein_ids ~ "Correlation",
    TRUE ~ "All"
  ))

# Calculate summary stats: median and SD
summary_stats <- expr_with_group %>%
  group_by(group) %>%
  summarise(
    median_expr = median(log2_mean_plus1, na.rm = TRUE),
    sd_expr = sd(log2_mean_plus1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Prepare the label text
summary_stats <- summary_stats %>%
  mutate(label = paste0("Median = ", round(median_expr, 2), "\nSD = ", round(sd_expr, 2)))

# Define custom colors
custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)

# Reorder factor levels
expr_with_group$group <- factor(expr_with_group$group, 
                                levels = c("All", "Correlation", "Unique"))
summary_stats$group <- factor(summary_stats$group, 
                              levels = c("All", "Correlation", "Unique"))

# Plot with all three groups
ggplot(expr_with_group, aes(x = group, y = log2_mean_plus1, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("All", "Correlation"), c("All", "Unique"), c("Correlation", "Unique")),
    label = "p.format",
  ) +
  geom_text(
    data = summary_stats,
    aes(x = group, y = max(expr_with_group$log2_mean_plus1) * 0.95, label = label),
    inherit.aes = FALSE,
    size = 4) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "none"
  )
