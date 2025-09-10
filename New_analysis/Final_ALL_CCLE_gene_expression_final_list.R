# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(cowplot)

# === 1. LOAD DATA ===
expr_df <- read.csv("Data/CCLE_expression.csv")
model_ID <-  read.csv("Data/cell_line_disease.csv")
# Rename the first column
colnames(expr_df)[1] <- "ModelID"
# Merge with model_ID data frame by "ModelID"
merged_df <- merge(model_ID, expr_df, by = "ModelID")
colnames(merged_df)[4:ncol(merged_df)] <- gsub("\\.\\..*?\\.", "", colnames(merged_df)[4:ncol(merged_df)])
gene_list <- read.csv("Data/final_gene_list.txt")

# Subset merged_df to keep only selected genes and metadata columns
filtered_df <- merged_df[, c("ModelID", "CellLineName", "DepmapModelType", intersect(gene_list$Gene, colnames(merged_df)))]
TCGA_modeltype <- read.csv("Data/TCGA_DepmapModelType.txt")


TCGA_modeltype_unique <- unique(trimws(as.character(TCGA_modeltype$DepmapModelType)))
# Clean up the column for matching
filtered_df$DepmapModelType <- trimws(as.character(filtered_df$DepmapModelType))

# Filter to keep only matching model types
filtered_df <- filtered_df[filtered_df$DepmapModelType %in% TCGA_modeltype_unique, ]




# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Create output directories
dir.create("GenePlots", showWarnings = FALSE)
dir.create("RankSummaries", showWarnings = FALSE)

# Color palette: SKCM = pink, UVM = purple, others = gray
depmap_colors <- function(types) {
  sapply(types, function(x) {
    if (x == "SKCM") "pink"
    else if (x == "UVM") "purple"
    else "gray"
  })
}

# Gather expression data into long format
long_df <- filtered_df %>%
  pivot_longer(cols = 4:ncol(filtered_df), names_to = "Gene", values_to = "Expression")

# Rank stats output
rank_stats_list <- list()

# Loop over unique genes
for (gene in unique(long_df$Gene)) {
  
  gene_df <- long_df %>%
    filter(Gene == gene) %>%
    mutate(DepmapModelType = fct_reorder(DepmapModelType, Expression, .fun = median, .desc = FALSE))  # sort by median
  
  # Plot
  p <- ggplot(gene_df, aes(x = DepmapModelType, y = Expression, fill = DepmapModelType)) +
    geom_boxplot(outlier.size = 0.8) +
    scale_fill_manual(values = depmap_colors(levels(gene_df$DepmapModelType))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Expression of", gene)) +
    ylab("Expression") + xlab("Model Type")
  
  # Save PDF
  ggsave(filename = paste0("GenePlots/", gene, "_boxplot.pdf"), plot = p, width = 8, height = 5)
  
  # Mean/Median summary
  stats <- gene_df %>%
    group_by(DepmapModelType) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      median_expr = median(Expression, na.rm = TRUE)
    ) %>%
    mutate(Gene = gene)
  
  rank_stats_list[[gene]] <- stats
}

# Combine and export stats
rank_summary_df <- bind_rows(rank_stats_list)
write.csv(rank_summary_df, "RankSummaries/gene_expression_summary.csv", row.names = FALSE)


# === Compute mean expression per gene per model type ===
ranked_df <- long_df %>%
  dplyr::group_by(Gene, DepmapModelType) %>%
  dplyr::summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(rank = rank(-mean_expr, ties.method = "average")) %>%
  dplyr::ungroup()

# === Summarize ranks across genes per model type ===
rank_summary <- ranked_df %>%
  dplyr::group_by(DepmapModelType) %>%
  dplyr::summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  )

# === Save rank summary ===
write.csv(rank_summary, "RankSummaries/gene_list_rank_summary.csv", row.names = FALSE)
saveRDS(rank_summary, "RankSummaries/gene_list_rank_summary.rds")

# Reshape for plotting
rank_long <- rank_summary %>%
  pivot_longer(cols = c(mean_rank, median_rank),
               names_to = "statistic",
               values_to = "rank_value") %>%
  dplyr::mutate(
    highlight = ifelse(DepmapModelType %in% c("SKCM", "UVM"), "highlight", "other"),
    fill_key = paste(statistic, highlight, sep = ".")
  )

# Define fill colors
color_map <- c(
  "mean_rank.highlight" = "#ff1493",    # SKCM/UVM mean rank - bright pink
  "median_rank.highlight" = "#8a2be2",  # SKCM/UVM median rank - bright purple
  "mean_rank.other" = "gray70",
  "median_rank.other" = "gray50"
)

# Create the plot
p_rank <- ggplot(rank_long, aes(x = reorder(DepmapModelType, rank_value), y = rank_value, fill = fill_key)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Mean and Median Rank of Gene List Expression by Model Type",
    x = "Model Type",
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

# Save plot
ggsave("CCLE_GENE_LIST_RANK_SUMMARY.PDF", plot = p_rank, width = 10, height = 6)



####
saveRDS(filtered_df,"CCLE_gene_disease_Figure1B.rds")

library(ggplot2)
library(dplyr)
library(forcats)  # for fct_reorder
## 
CCLE_gene_for_figure <- readRDS("CCLE_gene_disease_Figure1B.rds")
# Prepare the data
abr_df <- CCLE_gene_for_figure %>%
  dplyr::select(DepmapModelType = 3, ABR) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ABR, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ABR, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "blue", "Other" = "gray70")) +
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
