library(readr)
CCLE_trans_expr <- read_tsv("Data/CCLE_Melanoma_transcript_expression.txt")

transcript_ratio <- read.csv("Data/transcript_ratio.csv")
  
correlation <- read.csv("Data/NEW_correlation_protein_coding.csv")
final <- read.csv("Data/NEW_discordant_protein_coding.csv")

# Load libraries
library(dplyr)
library(ggplot2)

# Step 1: Clean transcript IDs in transcript_ratio
transcript_ratio <- transcript_ratio %>%
  mutate(transcript_id_clean = sub("\\..*$", "", transcript_id))

# Step 2: Clean correlation and final IDs (in case they aren't already clean)
correlation_ids <- sub("\\..*$", "", correlation$transcript_id)
final_ids <- sub("\\..*$", "", final$transcript_id)

# Step 3: Annotate group
transcript_ratio_grouped <- transcript_ratio %>%
  mutate(group = case_when(
    transcript_id_clean %in% final_ids ~ "Final",
    transcript_id_clean %in% correlation_ids ~ "Correlation",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Step 4: Boxplot only (no violin)
ggplot(transcript_ratio_grouped, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(
    title = "Transcript Ratio Comparison: Correlation vs Final",
    x = "Group",
    y = "Ratio"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(ratio ~ group, data = transcript_ratio_grouped)

# Format p-value nicely
p_val <- signif(wilcox_result$p.value, 3)
p_text <- paste0("Wilcoxon p = ", p_val)

# Plot with annotation
ggplot(transcript_ratio_grouped, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(
    title = "Transcript Ratio Comparison",
    x = "Group",
    y = "Ratio",
    subtitle = p_text
  ) +
  theme_minimal() +
  theme(legend.position = "none")

summary_stats <- transcript_ratio_grouped %>%
  group_by(group) %>%
  summarise(
    N = n(),
    Median = median(ratio, na.rm = TRUE),
    SD = sd(ratio, na.rm = TRUE)
  )
print(summary_stats)

# Prepare annotation text
summary_stats <- summary_stats %>%
  mutate(
    label = paste0("Median = ", round(Median, 2), 
                   "\nSD = ", round(SD, 2))
  )

# Create the plot
ggplot(transcript_ratio_grouped, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_text(
    data = summary_stats,
    aes(x = group, y = max(transcript_ratio_grouped$ratio) * 0.95, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 0
  ) +
  labs(
    title = "Transcript Ratio Comparison",
    x = "Group",
    y = "Ratio",
    subtitle = paste0("Wilcoxon p = ", signif(wilcox_result$p.value, 3))
  ) +
  theme_minimal() +
  theme(legend.position = "none")
