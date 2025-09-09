annotations <- readRDS("annotation_codes/hg38_annotated_transcripts_EXPANDING.rds")
final_list <- read.csv("Data/final_list_transcript_id.txt")

Final_list_annotated <- merge(final_list, annotations, by="transcript_id")




Final_gene_list <- read.csv("Data/final_gene_list.txt")
Gene_transcript_pearson <- read.csv("Data/merged_pearsons_annotated.csv")
# Replace NA Gene_name values with "Unknown"
Gene_transcript_pearson$Gene_name[is.na(Gene_transcript_pearson$Gene_name)] <- "Unknown"

# Create a highlight column
Gene_transcript_pearson$highlight <- ifelse(Gene_transcript_pearson$Gene_name == "MED15", "MED15", "Other"


# Remove everything from the first dot onwards in 'transcript_id' and 'gene_id'
Gene_transcript_pearson$transcript_id <- sub("\\..*$", "", Gene_transcript_pearson$transcript_id)
Gene_transcript_pearson$gene_id <- sub("\\..*$", "", Gene_transcript_pearson$gene_id)


library(ggplot2)

# Create a new column to flag whether the gene is MED15
Gene_transcript_pearson$highlight <- ifelse(Gene_transcript_pearson$Gene_name == "STARD10", "STARD10", "Other")

# Separate MED15 and other points
med15_df <- subset(Gene_transcript_pearson, Gene_name == "STARD10")
other_df <- subset(Gene_transcript_pearson, Gene_name != "STARD10")

# Plot: plot 'other' first, then 'MED15' on top
ggplot() +
  geom_point(data = other_df, aes(x = pearson_transcript, y = pearson_gene),
             color = "gray30", size = 0.025, alpha = 0.3) +
  geom_point(data = med15_df, aes(x = pearson_transcript, y = pearson_gene),
             color = "red", size = 2, alpha = 0.9) +
  labs(
    title = "Pearson Correlation: Transcript vs Gene",
    x = "Transcript Pearson",
    y = "Gene Pearson"
  ) +
  theme_minimal()
  
  
  ### Plot for all genes
  # Start PDF output
pdf("CCLE_Pearson_Gene_transcript.pdf", width = 7, height = 6)

# Loop through each gene in Final_gene_list
for (gene in Final_gene_list$Gene_name) {
  
  # Subset data
  highlight_df <- subset(Gene_transcript_pearson, Gene_name == gene)
  other_df <- subset(Gene_transcript_pearson, Gene_name != gene)
  
  # Create plot
  p <- ggplot() +
    geom_point(data = other_df, aes(x = pearson_transcript, y = pearson_gene),
               color = "gray30", size = 0.025, alpha = 0.3) +
    geom_point(data = highlight_df, aes(x = pearson_transcript, y = pearson_gene),
               color = "red", size = 2, alpha = 0.9) +
    labs(
      title = paste("Pearson Correlation: Transcript vs Gene\nHighlight:", gene),
      x = "Transcript Pearson",
      y = "Gene Pearson"
    ) +
    theme_minimal()
  
  # Print plot to PDF
  print(p)
}

# Close the PDF device
dev.off()



### PLOT FOR ABR FOR PAPER
library(ggplot2)
library(dplyr)
# Annotate the points
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    Gene_name == "ABR" & transcript_id == "ENST00000544583" ~ "ABR_Main",
    Gene_name == "ABR" ~ "ABR_Other",
    TRUE ~ "Other"
  ))

library(ggplot2)
library(dplyr)

# Annotate the points
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    Gene_name == "ABR" & transcript_id == "ENST00000544583" ~ "ABR_Main",
    Gene_name == "ABR" ~ "ABR_Other",
    TRUE ~ "Other"
  ))

# Plot with updated column names
ggplot() +
  # Background points (gray, small, transparent)
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Other"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "gray30", size = 0.1, alpha = 0.3
  ) +
  # Highlighted ABR transcripts
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight != "Other"),
    aes(x = pearson_gene, y = pearson_transcript, color = highlight),
    size = 2.5
  ) +
  scale_color_manual(values = c("ABR_Main" = "blue", "ABR_Other" = "red")) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()
  
  
  
  library(ggplot2)
library(dplyr)
Gene_transcript_pearson <- read.csv("Data/merged_pearsons_annotated.csv")
Gene_transcript_pearson$transcript_id <- sub("\\..*$", "", Gene_transcript_pearson$transcript_id)
Gene_transcript_pearson$gene_id <- sub("\\..*$", "", Gene_transcript_pearson$gene_id)



# Define the transcript IDs
blue_transcript <- "ENST00000544583"

red_transcripts <- c(
  "ENST00000291107", "ENST00000302538", "ENST00000536794", "ENST00000543210",
  "ENST00000570441", "ENST00000570525", "ENST00000570688", "ENST00000571022",
  "ENST00000571120", "ENST00000571306", "ENST00000571383", "ENST00000571543",
  "ENST00000571797", "ENST00000571945", "ENST00000572152", "ENST00000572441",
  "ENST00000572585", "ENST00000572650", "ENST00000573325", "ENST00000573559",
  "ENST00000573667", "ENST00000573895", "ENST00000574048", "ENST00000574139",
  "ENST00000574257", "ENST00000574266", "ENST00000574437", "ENST00000574544",
  "ENST00000574632", "ENST00000574875", "ENST00000575770", "ENST00000575934",
  "ENST00000576668", "ENST00000576964", "ENST00000577052"
)

# Add highlight column
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    transcript_id == blue_transcript ~ "Blue",
    transcript_id %in% red_transcripts ~ "Red",
    TRUE ~ "Other"
  ))

# Plot
ggplot() +
  # Background points (small and transparent)
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Other"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "gray70", size = 0.1, alpha = 0.3
  ) +
  # Highlighted points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight != "Other"),
    aes(x = pearson_gene, y = pearson_transcript, color = highlight),
    size = 2.0
  ) +
  scale_color_manual(values = c("Blue" = "red", "Red" = "blue")) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()

