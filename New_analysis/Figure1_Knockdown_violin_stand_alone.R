transcript_annotation <- readRDS("hg38_annotated_transcripts_EXPANDING.rds")

transcript_unique <- read.csv("Data/Final_CCLE_TCGA_overlap.csv")
transcript_correlation_all <-read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")


library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Step 1: Clean knockdown data
knockdown_data <- transcript_annotation %>%
  filter(!is.nan(avg_knockdown)) %>%
  filter(avg_knockdown > 0 & avg_knockdown <= 2) %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  select(transcript_id, group, avg_knockdown)

# Custom group fill colors
custom_colors <- c(
  "All" = "#FADADD",
  "Correlation" = "#F08080",
  "Unique" = "#FF0000"
)

# Step 2: Run statistical tests
kruskal_test <- compare_means(avg_knockdown ~ group, data = knockdown_data, method = "kruskal.test")

pairwise_tests <- compare_means(
  avg_knockdown ~ group,
  data = knockdown_data,
  method = "wilcox.test",
  p.adjust.method = "fdr"
)

# Add feature label for later merging
kruskal_test$Feature <- "avg_knockdown"
pairwise_tests$Feature <- "avg_knockdown"

# Step 3: Plot violin + boxplot with stat annotations
ggplot(knockdown_data, aes(x = group, y = avg_knockdown, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(method = "kruskal.test", label.y = 2.1) +  # Global
  stat_compare_means(
    comparisons = list(
      c("All", "Correlation"),
      c("All", "Unique"),
      c("Correlation", "Unique")
    ),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = c(1.9, 2.0, 2.1)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Filtered avg_knockdown per Transcript Group (0 < x ≤ 2)",
    y = "Average Knockdown"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, margin = margin(r = 15))
  )

# Step 4: Format and export the result table
pairwise_table <- pairwise_tests %>%
  select(Feature, group1, group2, method, p, p.adj, p.signif)

# Save to CSV
write.csv(pairwise_table, "avg_knockdown_pairwise_stats.csv", row.names = FALSE)
