correlations <- read.csv("Data/correlation_list.csv")


# Step 1: Subset individual correlation columns
tcga_spearman_hits <- correlations[!is.na(correlations$spearman_TCGA_value) & correlations$spearman_TCGA_value > 0.3, "transcriptid_spearman_TCGA"]
tcga_pearson_hits  <- correlations[!is.na(correlations$pearson_TCGA_value)  & correlations$pearson_TCGA_value  > 0.3, "transcriptid_pearson_TCGA"]

ccle_spearman_hits <- correlations[!is.na(correlations$spearman_CCLE_value) & correlations$spearman_CCLE_value > 0.5, "transcriptid_spearman_CCLE"]
ccle_pearson_hits  <- correlations[!is.na(correlations$pearson_CCLE_value)  & correlations$pearson_CCLE_value  > 0.5, "transcriptid_pearson_CCLE"]

# Step 2: Combine lists within datasets (union)
tcga_transcripts <- unique(c(tcga_spearman_hits, tcga_pearson_hits))
ccle_transcripts <- unique(c(ccle_spearman_hits, ccle_pearson_hits))

# Step 3: Final intersection (only transcripts present in both)
final_transcripts_new <- intersect(tcga_transcripts, ccle_transcripts)

# Step 4: Write to file
write.csv(data.frame(transcript_id = final_transcripts_new), "correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv", row.names = FALSE)
