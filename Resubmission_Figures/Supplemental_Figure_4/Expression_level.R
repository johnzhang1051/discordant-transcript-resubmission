library(readr)
library(dplyr)
library(ggplot2)

CCLE_trans_expr <- read_tsv("data/CCLE_Melanoma_transcript_expression.txt")
CCLE_trans_expr <- CCLE_trans_expr %>%
  mutate(transcript_id = sub("\\..*$", "", transcript_id))
CCLE_trans_expr_clean <- CCLE_trans_expr[,-1]


# Step 1: Calculate the row means (average counts per transcript)
CCLE_trans_expr_clean$mean_count <- rowMeans(CCLE_trans_expr_clean[, 2:64], na.rm = TRUE)

# Step 2: Log2 transform (mean + 1)
CCLE_trans_expr_clean$log2_mean_plus1 <- log2(CCLE_trans_expr_clean$mean_count + 1)

CCLE_trans_expr_log2 <- CCLE_trans_expr_clean[,-(2:65)]


correlation <- read.csv("data/correlated_RESUBMISSION.csv")
final <- read.csv("data/discordant_RESUBMISSION.csv")

## plot average gene expression comparison
# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Assume CCLE_trans_expr_log2, correlation, and final already loaded
# Make sure 'correlation' and 'final' are vectors of transcript_ids
correlation_ids <- unique(correlation$transcript_id)
final_ids <- unique(final$transcript_id)

# Assign group labels
expr_with_group <- CCLE_trans_expr_log2 %>%
  mutate(group = case_when(
    transcript_id %in% final_ids ~ "Discordant",
    transcript_id %in% correlation_ids ~ "Correlation",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Calculate summary stats: median and SD
summary_stats <- expr_with_group %>%
  group_by(group) %>%
  summarise(
    median_expr = median(log2_mean_plus1, na.rm = TRUE),
    sd_expr = sd(log2_mean_plus1, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare the label text (format: "Median ± SD")
summary_stats <- summary_stats %>%
  mutate(label = paste0("Median = ", round(median_expr, 2), "\nSD = ", round(sd_expr, 2)))

# Plot


p <- ggplot(expr_with_group, aes(x = group, y = log2_mean_plus1, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    comparisons = list(c("Correlation", "Discordant")),
    label.y = max(expr_with_group$log2_mean_plus1) * 1.05  # slightly above top
  ) +
  geom_text(
    data = summary_stats,
    aes(x = group, y = max(expr_with_group$log2_mean_plus1) * 0.95, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = "Comparison of Average Transcript Expression",
    x = "Group",
    y = "Mean log2(Expression + 1)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Show the plot
print(p)

custom_colors <- c(
  "Correlation" = "#F08080",  # medium pink
  "Discordant" = "#FF0000"         # red
)

p <- ggplot(expr_with_group, aes(x = group, y = log2_mean_plus1, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    comparisons = list(c("Correlation", "Discordant")),
    label.y = max(expr_with_group$log2_mean_plus1, na.rm = TRUE) * 1.05
  ) +
  geom_text(
    data = summary_stats,
    aes(x = group, y = max(expr_with_group$log2_mean_plus1, na.rm = TRUE) * 0.95, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Comparison of Average Transcript Expression",
    x = "Group",
    y = "Mean log2(Expression + 1)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

