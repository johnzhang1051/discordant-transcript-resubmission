final_list <- read.csv("Data/Final_CCLE_TCGA_overlap.csv")
correlations_list  <- read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")
gene_transcript <-  read.csv("Data/gene_transcript.csv")


library(dplyr)

# Merge final_list
final_list_annotated <- final_list %>%
  left_join(gene_transcript, by = "transcript_id")

# Merge correlations_list
correlations_list_annotated <- correlations_list %>%
  left_join(gene_transcript, by = "transcript_id")

write.csv(final_list_annotated, "final_list_new.csv")
write.csv(correlations_list_annotated, "correlations_list_new.csv")



# Load required packages
library(org.Hs.eg.db)
library(AnnotationDbi)

final_list_genes <- read.csv("final_list_genes_new.csv")
correlation_list_genes <- read.csv("correlations_list_genes_new.csv")

# Map Gene Symbols to Ensembl Gene IDs
correlation_list_genes$gene_id <- mapIds(
  org.Hs.eg.db,
  keys = correlation_list_genes$Gene,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"  # or use "list" if you want all matches
)


write.csv(final_list_genes,"final_list_genes_geneID.csv")
write.csv(correlation_list_genes, "correlation_list_genes_geneID.csv")