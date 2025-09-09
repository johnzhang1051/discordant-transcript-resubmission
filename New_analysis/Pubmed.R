# Load required packages
install.packages("easyPubMed")
install.packages("tidyverse")

library(easyPubMed)
library(tidyverse)

# Load gene lists
final_list <- read.csv("Data/final_list_genes_new.csv")
correlation_list <- read.csv("Data/correlations_list_genes_new.csv")

# Extract gene names
genes_final <- unique(final_list$Gene_name)
genes_correlation <- unique(correlation_list$Gene_name)

# Function to search PubMed for MITF-related mentions
get_pubmed_hits <- function(gene_symbol) {
  query <- paste0(gene_symbol, " AND MITF AND (regulation OR binding OR target)")
  result <- get_pubmed_ids(query)
  return(result$Count)
}

# Apply function to each gene list
hits_final <- data.frame(
  Gene = genes_final,
  Hits = sapply(genes_final, get_pubmed_hits),
  Group = "Final"
)

hits_correlation <- data.frame(
  Gene = genes_correlation,
  Hits = sapply(genes_correlation, get_pubmed_hits),
  Group = "Correlation"
)

# Combine into one dataframe
# Combine into one dataframe
# Combine and ensure Hits is numeric
all_hits <- bind_rows(hits_final, hits_correlation) %>%
  mutate(Hits = as.numeric(Hits))

# Summary stats
summary_stats <- all_hits %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Hits),
    Median = median(Hits),
    .groups = "drop"
  )


library(ggplot2)
library(dplyr)
library(scales)

# Set custom colors for the two groups
custom_colors <- c("Final" = "#D73027", "Correlation" = "#F08080")  # red and medium pink

# Plot
all_hits %>%
  mutate(has_pubmed_hit = Hits > 1) %>%
  group_by(Group) %>%
  summarise(prop_hit = mean(has_pubmed_hit), .groups = "drop") %>%
  ggplot(aes(x = Group, y = prop_hit, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.85) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.15)
  ) +
  labs(
    title = "Proportion of Genes with >1 MITF-Related PubMed Hit",
    y = "Percent with PubMed Evidence",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
