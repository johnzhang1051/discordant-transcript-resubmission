library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Load data
CCLE_trans_expr <- read_tsv("data/CCLE_Melanoma_transcript_expression.txt")
CCLE_trans_expr <- CCLE_trans_expr %>%
  mutate(transcript_id = sub("\\..*$", "", transcript_id))

CCLE_trans_expr_clean <- CCLE_trans_expr[,-1]

#### FILTER TO ONLY MITF-HIGH EXPRESSING CELL LINES:
# This means only cell-lines with MITF expression > median expression of MITF
# Get MITF transcript expression (ENST00000394351)
mitf_expr <- CCLE_trans_expr_clean %>%
  filter(transcript_id == "ENST00000394351")

# Extract the expression values as a numeric vector
mitf_values <- as.numeric(mitf_expr[1, -1])

# Reshape to long format for plotting
mitf_long <- mitf_expr %>%
  select(-transcript_id) %>%
  pivot_longer(cols = everything(), 
               names_to = "cell_line", 
               values_to = "expression")

# Calculate median for threshold line
mitf_median <- median(mitf_long$expression, na.rm = TRUE)

# Add classification
mitf_long <- mitf_long %>%
  mutate(classification = ifelse(expression >= mitf_median, "MITF-High", "MITF-Low"))

# Create histogram
ggplot(mitf_long, aes(x = expression, fill = classification)) +
  geom_histogram(bins = 20, alpha = 0.7, color = "black") +
  geom_vline(xintercept = mitf_median, 
             color = "red", linetype = "dashed", size = 1.2) +
  scale_fill_manual(values = c("MITF-High" = "darkred", "MITF-Low" = "steelblue")) +
  annotate("text", x = mitf_median, y = Inf, 
           label = paste("Median =", round(mitf_median, 2)), 
           color = "red", vjust = 2, hjust = -0.1, fontface = "bold") +
  labs(
    title = "MITF Expression Distribution in Melanoma Cell Lines",
    subtitle = "Transcript: ENST00000394351",
    x = "MITF Expression (TPM)",
    y = "Number of Cell Lines",
    fill = "Classification"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "bottom"
  )

# Filter to only MITF-high cell lines
mitf_high_cells <- names(mitf_expr[1, -1])[mitf_values >= mitf_median]

CCLE_trans_expr_clean <- CCLE_trans_expr_clean %>%
  select(transcript_id, all_of(mitf_high_cells))

########## NOW WE RUN THE ANALYSIS AS USUAL ##########


# Step 1: Calculate the row means (average counts per transcript)
CCLE_trans_expr_clean$mean_count <- rowMeans(CCLE_trans_expr_clean[, -1], na.rm = TRUE)

# Step 2: Log2 transform (mean + 1)
CCLE_trans_expr_clean$log2_mean_plus1 <- log2(CCLE_trans_expr_clean$mean_count + 1)

CCLE_trans_expr_log2 <- CCLE_trans_expr_clean %>%
  select(transcript_id, log2_mean_plus1)

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
