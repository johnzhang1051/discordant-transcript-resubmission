# If not already installed
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

####### Get Data:

# Pull discordant transcripts
final_protein_coding <- read.csv("Data/discordant_RESUBMISSION.csv")

# Pull correlated transcripts
correlatin_coding_only <- read.csv("Data/correlated_RESUBMISSION.csv")


# Load transcript_to_gene mapping
transcript_gene_map <- read.csv("Data/transcript_to_gene.csv")
transcript_gene_map$transcript_id <- sub("\\..*$", "", transcript_gene_map$transcript_id)

final_protein_GENE <- merge(final_protein_coding, transcript_gene_map, by='transcript_id')
#final_protein_GENE <- as.data.frame(final_protein_GENE[,-(1:3)])

correlation_protein_GENE <- merge(correlatin_coding_only, transcript_gene_map, by='transcript_id')
#correlation_protein_GENE <- as.data.frame(correlation_protein_GENE[,-(1:3)])

write.csv(final_protein_GENE,"Data/final_protein_coding_GENES.csv", row.names = FALSE)
write.csv(correlation_protein_GENE,"Data/correlation_protein_coding_GENES.csv", row.names = FALSE)


#all_protein_coding <- read.csv("Data/protein_coding_RESUBMISSION.csv")
#all_protein_GENE <- merge(all_protein_coding, transcript_gene_map, by='transcript_id')

###################### Start Graphing
correlation_genes <- read.csv("Data/correlation_protein_coding_GENES.csv")
names(correlation_genes)[names(correlation_genes) == "Gene.y"] <- "Gene_name"


final_genes <- read.csv ("Data/final_protein_coding_GENES.csv")
names(final_genes)[names(final_genes) == "Gene.y"] <- "Gene_name"
names(final_genes)[names(final_genes) == "Gene"] <- "Gene_name"
names(final_genes)[names(final_genes) == "Gene.name"] <- "Gene_name"

## extract to vector
correlation_symbols <- as.character(correlation_genes$Gene_name)
final_symbols <- as.character(final_genes$Gene_name)

# Now run bitr with the vectors
correlation_entrez <- bitr(correlation_symbols, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)

final_entrez <- bitr(final_symbols, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

## GO terms
# Correlation Genes
go_corr <- enrichGO(gene         = correlation_entrez$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.05,
                    readable     = TRUE)

# Final Genes
go_final <- enrichGO(gene         = final_entrez$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod= "BH",
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     readable     = TRUE)

# View top results
head(go_corr)
head(go_final)

library(ggplot2)
dotplot(go_corr, showCategory = 15) +
  ggtitle("GO Biological Processes Enriched in correlation_final")
barplot(go_corr, showCategory = 15, title = "Top GO Terms")

dotplot(go_final, showCategory = 15) +
  ggtitle("GO Biological Processes Enriched in discordant_final")
barplot(go_final, showCategory = 15, title = "Top GO Terms")


# Save to CSV
write.csv(as.data.frame(go_corr), "go_results/GO_correlation_genes.csv", row.names = FALSE)
write.csv(as.data.frame(go_final), "go_results/GO_final_genes.csv", row.names = FALSE)

## genes linked to specific GO
library(dplyr)
library(stringr)
# Search for the GO term (flexible search)
go_of_interest <- go_final %>%
  filter(str_detect(Description, "positive regulation of receptor clustering"))

# See genes linked
genes_linked <- strsplit(go_of_interest$geneID[1], "/")[[1]]

# View them
genes_linked






# Reactome for correlation_genes
reactome_correlation <- enrichPathway(gene = correlation_entrez$ENTREZID,
                                pvalueCutoff = 0.05,
                                readable = TRUE)
head(reactome_correlation)


# Reactome for final_genes
reactome_final <- enrichPathway(gene = final_entrez$ENTREZID,
                                pvalueCutoff = 0.05,
                                readable = TRUE)

# View Reactome results
head(reactome_final)
## ignore signficiance
reactome_final_all <- enrichPathway(
  gene         = final_entrez$ENTREZID,
  organism     = "human",
  pvalueCutoff = 0.5, # include all results
  qvalueCutoff = 0.5, # include all results
  readable     = TRUE
)
# Convert to data frame
reactome_df <- as.data.frame(reactome_final_all)

# View top 10 by raw p-value
head(reactome_df[order(reactome_df$pvalue), ], 10)

# Or view top 10 by z-score (if you want ranked activity)
head(reactome_df[order(-reactome_df$zScore), ], 10)
##Even if they're not significant, you can still visualize them:
dotplot(reactome_final_all, showCategory = 15) +
  ggtitle("Top Reactome Pathways (Unfiltered)")


library(ReactomePA)
dotplot(reactome_correlation, showCategory = 15) +
  ggtitle("Reactome Pathways Enriched in correlation_final")

# Add source label
## update if you want to merge lists
go_df <- as.data.frame(go_result) %>% mutate(Source = "GO")
reactome_df <- as.data.frame(reactome_result) %>% mutate(Source = "Reactome")
# Combine
combined_enrichment <- bind_rows(go_df, reactome_df)
# Save combined results
write.csv(combined_enrichment, "go_results.Combined_GO_Reactome_correlation_final.csv", row.names = FALSE)




## KEGG Pathway
# KEGG requires organism code — hsa = human
kegg_corr <- enrichKEGG(gene = correlation_entrez$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)

kegg_final <- enrichKEGG(gene = final_entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.5)

# Dotplot for KEGG results
kegg_corr@result$Group <- "Correlation Genes"
kegg_final@result$Group <- "Final Genes"
combined_KEGG <- merge_result(list(correlation = kegg_corr, 
                                   discordant = kegg_final))
dotplot(combined_KEGG, showCategory = 10, split = "Group") +
  ggtitle("KEGG Pathway Enrichment")



### how to tell which reactome hits 
reactome_correlation@result$Description

# Filter for "MITF-M regulated melanocyte development"
mitf_pathway <- reactome_correlation@result %>%
  dplyr::filter(Description == "Ion channel transport")

# Extract the gene IDs (they are usually concatenated with / or ;)
genes_in_mitf_pathway <- strsplit(mitf_pathway$geneID, "/")[[1]]  # If / separated
# If semicolon-separated, use ";" instead

# View the genes
genes_in_mitf_pathway
