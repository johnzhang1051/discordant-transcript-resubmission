library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

expr_df <- readRDS("Data/CCLE_Melanoma_transcript_expression.rds")
# Remove trailing period from transcript_id
expr_df$transcript_id <- sub("\\.$", "", expr_df$transcript_id)
# Step 1: Read in final_list.txt from the 'Data' folder
final_list <- read.csv("Data/final_list.csv")
# Step 5: Merge expr_df with final_list by 'gene_id'
merged_df <- merge(expr_df, final_list, by = "gene_id")
colnames(merged_df)[colnames(merged_df) == "transcript_id.x"] <- "transcript_id"

# Assuming:
# merged_df$transcript_id is the transcript ID column
# final_list$V1 contains the list of transcript IDs to match

# Subset for matching transcripts
transcript_unique <- merged_df[merged_df$transcript_id %in% final_list[[1]], ]
# Copy the data frame to avoid modifying original
transcript_unique_log <- transcript_unique
# Apply log2(x + 1) to expression values only
transcript_unique_log[, 3:65] <- log2(transcript_unique_log[, 3:65] + 1)


# Subset for non-matching transcripts
transcript_all_other <- merged_df[!merged_df$transcript_id %in% final_list[[1]], ]
library(dplyr)

# Identify expression columns by excluding transcript_id and gene_id
expr_cols <- names(transcript_all_other)[!(names(transcript_all_other) %in% c("transcript_id", "gene_id"))]

# Group by gene_id and sum numeric expression values
transcript_all_other_summed <- transcript_all_other %>%
  group_by(gene_id) %>%
  summarise(across(all_of(expr_cols), ~ sum(as.numeric(.), na.rm = TRUE)), .groups = 'drop')

# Log2(x + 1) transform
transcript_all_other_log <- transcript_all_other_summed
transcript_all_other_log[, expr_cols] <- log2(transcript_all_other_log[, expr_cols] + 1)



mitf_exp <- readRDS("Data/MITF_transcript.rds")
mitf_exp <- as.data.frame(mitf_exp)
colnames(mitf_exp)[colnames(mitf_exp) == "ENST00000394351."] <- "mitf"

# put in same order
# Step 1: Extract sample names from expr_df (columns 3 to 65)
expr_sample_names <- colnames(expr_df)[3:65]

# Step 2: Reorder mitf_exp rows to match expr_df sample order
mitf_exp_ordered <- mitf_exp[expr_sample_names, , drop = FALSE]

# Step 3: (Optional) Check alignment
all(rownames(mitf_exp_ordered) == expr_sample_names)

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Prepare sample ID map
sample_ids <- colnames(transcript_unique_log)[3:65]
mitf_df <- data.frame(sample = sample_ids, mitf = mitf_exp)
mitf_df$sample <- as.character(mitf_df$sample)

# Store correlation results
correlation_results <- data.frame()

# Loop over genes
unique_genes <- unique(final_list$gene_id)

for (gene in unique_genes) {
  
  # --- Get transcript_unique expression ---
  tu_expr <- transcript_unique_log %>%
    filter(gene_id == gene) %>%
    select(3:65) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "gene_expr") %>%
    mutate(group = "transcript_unique")
  
  # --- Get transcript_all_other expression ---
  tao_expr <- transcript_all_other_log %>%
    filter(gene_id == gene) %>%
    pivot_longer(cols = 2:64, names_to = "sample", values_to = "gene_expr") %>%
    mutate(group = "transcript_all_other")
  
  # Combine both
  combined_expr <- bind_rows(tu_expr, tao_expr)
  combined_expr$sample <- as.character(combined_expr$sample)
  
  # Merge with MITF expression
  combined_expr <- left_join(combined_expr, mitf_df, by = "sample")
  
  # Skip if something went wrong with merge
  if (!"mitf" %in% colnames(combined_expr)) next
  
  # --- Compute correlations by group ---
  cor_data <- combined_expr %>%
    group_by(group) %>%
    summarise(
      pearson = cor(mitf, gene_expr, method = "pearson"),
      spearman = cor(mitf, gene_expr, method = "spearman"),
      .groups = "drop"
    ) %>%
    mutate(gene_id = gene)
  
  # Save to cumulative results
  correlation_results <- bind_rows(correlation_results, cor_data)
  
  # --- Make label text for plot ---
  label_text <- cor_data %>%
    mutate(label = paste0(group, ": r = ", round(pearson, 2), 
                          ", ρ = ", round(spearman, 2))) %>%
    pull(label) %>%
    paste(collapse = "\n")
  
  
  # --- Plot ---
  p <- ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, label = label_text, hjust = -0.05, vjust = 1.2, size = 4, color = "black") +
    labs(title = paste("MITF vs", gene, "Expression"),
         x = "MITF expression (log2 + 1)",
         y = paste(gene, "expression (log2 + 1)")) +
    theme_minimal()
  
   print(p)
}

# --- Write correlations to file ---
write.csv(correlation_results, "MITF_gene_correlations.csv", row.names = FALSE)


#### PLOTTING ABR ONLY
library(dplyr)
library(tidyr)
library(ggplot2)

# Extract ABR expression data
gene <- "ABR"

# Expression for transcript_unique
tu_expr <- transcript_unique_log %>%
  filter(gene_id == gene) %>%
  select(3:65) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "gene_expr") %>%
  mutate(group = "transcript_unique")

# Expression for transcript_all_other
tao_expr <- transcript_all_other_log %>%
  filter(gene_id == gene) %>%
  pivot_longer(cols = 2:64, names_to = "sample", values_to = "gene_expr") %>%
  mutate(group = "transcript_all_other")

# Combine and merge MITF
combined_expr <- bind_rows(tu_expr, tao_expr)
combined_expr$sample <- as.character(combined_expr$sample)
combined_expr <- left_join(combined_expr, mitf_df, by = "sample")

# Compute correlations
cor_data <- combined_expr %>%
  group_by(group) %>%
  summarise(
    pearson = cor(mitf, gene_expr, method = "pearson"),
    spearman = cor(mitf, gene_expr, method = "spearman"),
    .groups = "drop"
  )

# Prepare label
label_text <- cor_data %>%
  mutate(label = paste0(group, ": r = ", round(pearson, 2),
                        ", ρ = ", round(spearman, 2))) %>%
  pull(label) %>%
  paste(collapse = "\n")

# Plot
ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  annotate("text", x = -Inf, y = Inf, label = label_text, hjust = -0.05, vjust = 1.2, size = 4, color = "black") +
  labs(title = paste("MITF vs", gene, "Expression"),
       x = "MITF expression (log2 + 1)",
       y = paste(gene, "expression (log2 + 1)")) +
  theme_minimal()


##
# Modern clean plot with better styling (ABR vs MITF)
library(ggplot2)

ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
  geom_point(size = 2, alpha = 0.9, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.2) +
  labs(
    title = "ABR Expression vs MITF",
    subtitle = "Colored by transcript group",
    x = "MITF expression (log2 + 1)",
    y = "ABR expression (log2 + 1)",
    color = "Transcript group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

### for paper plot
ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.1) +
  scale_color_manual(
    values = c("transcript_unique" = "blue", "transcript_all_other" = "red")
  ) +
  labs(
    title = "ABR vs MITF Expression",
    x = "MITF expression (log2 + 1)",
    y = "ABR expression (log2 + 1)",
    color = "Transcript group"
  ) +
  theme_classic(base_size = 13, base_family = "Times") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


# Create the plot and assign it to an object
abr_plot <- ggplot(combined_expr, aes(x = mitf, y = gene_expr, color = group)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.1) +
  scale_color_manual(
    values = c("transcript_unique" = "red", "transcript_all_other" = "blue")
  ) +
  labs(
    title = "ABR vs MITF Expression",
    x = "MITF expression (log2 + 1)",
    y = "ABR expression (log2 + 1)",
    color = "Transcript group"
  ) +
  theme_classic(base_size = 13, base_family = "Times") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save to PDF
ggsave("ABR_CCLE_GENE_TRANSCRIPT.PDF", plot = abr_plot, width = 7, height = 5)
