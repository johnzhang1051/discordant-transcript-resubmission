library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)


expr_df <- fread("Data/Transcript_expression_melanoma_log2.csv")
transcript_type <- read.csv("Data/transcripttype.csv")
# Ensure required package is loaded
library(dplyr)

# Step 1: Get vector of protein-coding transcript IDs
protein_coding_ids <- transcript_type %>%
  filter(transcript_type == "protein_coding") %>%
  pull(transcript_ID)

# Step 2: Subset expr_df
# Keep Sample_ID column and only those ENST columns present in protein_coding_ids
expr_df <- expr_df %>%
  select(Sample_ID, all_of(intersect(colnames(expr_df)[-1], protein_coding_ids)))






final_discordant <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")
transcript_to_gene <- read.csv("Data/human_transcript_gene_map.csv")
mitf_test<- expr_df %>%
  select(Sample_ID, ENST00000394351)

library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Clean up expression matrix
expr_mat <- expr_df
expr_mat$Sample_ID <- make.unique(as.character(expr_mat$Sample_ID))
expr_mat <- expr_mat %>% column_to_rownames("Sample_ID")

# Step 2: Filter transcript_to_gene to relevant transcripts
transcript_to_gene_filtered <- transcript_to_gene %>%
  filter(transcript_ID %in% colnames(expr_mat))

# Step 3: MITF expression vector (still log2(TPM+1))
mitf_id <- "ENST00000394351"
if (!mitf_id %in% colnames(expr_mat)) stop("MITF transcript not found in expr_df.")
mitf_expr <- expr_mat[[mitf_id]]

# Step 4: Back-transform to TPM
expr_tpm <- 2^expr_mat - 1

# Output folder
dir.create("MITF_transcript_vs_geneSummed", showWarnings = FALSE)

# Step 5: Loop over focal transcripts
for (i in seq_len(nrow(final_discordant))) {
  focal_tx <- final_discordant$transcript_ID[i]
  gene <- final_discordant$Gene[i]
  
  if (!focal_tx %in% colnames(expr_mat)) next
  
  # Get other transcripts from same gene
  all_gene_tx <- transcript_to_gene_filtered %>%
    filter(Gene == gene) %>%
    pull(transcript_ID)
  
  other_tx <- setdiff(all_gene_tx, focal_tx)
  valid_tx <- other_tx[other_tx %in% colnames(expr_tpm)]
  if (length(valid_tx) == 0) next
  
  # Extract expressions
  focal_expr <- expr_mat[[focal_tx]]  # already log2(TPM+1)
  summed_tpm <- rowSums(expr_tpm[, valid_tx, drop = FALSE], na.rm = TRUE)
  summed_expr <- log2(summed_tpm + 1)  # summed TPM → log2(TPM+1)
  
  # Prepare data frame for plotting
  plot_df <- data.frame(
    Sample_ID = rownames(expr_mat),
    MITF = mitf_expr,
    Focal_Transcript = focal_expr,
    Summed_Other_Transcripts = summed_expr
  ) %>%
    pivot_longer(cols = c("Focal_Transcript", "Summed_Other_Transcripts"),
                 names_to = "group", values_to = "gene_expr") %>%
    mutate(group = recode(group,
                          "Focal_Transcript" = "transcript_unique",
                          "Summed_Other_Transcripts" = "transcript_all_other"))
  
  # Plot in ABR-style
  abr_plot <- ggplot(plot_df, aes(x = MITF, y = gene_expr, color = group)) +
    geom_point(size = 2, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.1) +
    scale_color_manual(
      values = c("transcript_unique" = "red", "transcript_all_other" = "blue")
    ) +
    labs(
      title = paste(gene, "vs MITF Expression"),
      subtitle = paste("Focal transcript:", focal_tx),
      x = "MITF expression (log2 + 1)",
      y = paste(gene, "expression (log2 + 1)"),
      color = "Transcript group"
    ) +
    theme_classic(base_size = 13, base_family = "Times") +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # Display in RStudio
  print(abr_plot)
  
  # Save to PDF
  ggsave(
    filename = paste0("MITF_transcript_vs_geneSummed/", focal_tx, "_MITF_vs_geneSum.pdf"),
    plot = abr_plot, width = 7, height = 5
  )
}
