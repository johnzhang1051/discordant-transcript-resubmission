library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

expr_df <- readRDS("Data/Supplemental File 3 - TCGA_gene_expression.rds")
expr_df$sample <- sub("\\..*$", "", expr_df$sample)
# Move 'sample' column to row names
rownames(expr_df) <- expr_df$sample

# Remove the 'sample' column
expr_df$sample <- NULL


