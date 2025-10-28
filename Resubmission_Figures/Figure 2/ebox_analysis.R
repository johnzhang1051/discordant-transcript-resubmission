# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load your new data sources (instead of the old lists)
protein_coding <- read.csv("data/protein_coding_RESUBMISSION.csv")
transcript_unique <- read.csv("data/discordant_RESUBMISSION.csv")  # discordant transcripts
transcript_correlation_all <- read.csv("data/correlated_RESUBMISSION.csv")  # correlated transcripts

# Create GenomicRanges objects (same as original)
Eboxes_GR <- readGeneric("data/EBoxes.txt", chr = 1, start = 2, end = 3, strand = NULL,
                         meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                         remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

# Load transcript database and get promoter regions (same as original)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
transcripts <- transcripts(txdb)
prom <- promoters(transcripts, upstream = 2000, downstream = 200)

# Find EBOX overlaps with promoters (same as original)
EBOXProm.hits <- findOverlaps(Eboxes_GR, prom)
EBOX_prom <- transcripts[subjectHits(EBOXProm.hits)]
write.table(EBOX_prom, file = "EBOXPROM.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Count EBOX motifs per transcript (same as original)
EBOX_prom_2 <- read.csv(file = "EBOXPROM.csv", row.names = NULL, head = TRUE, sep = ",", stringsAsFactors = FALSE)
EBOX_transcripts_counts <- EBOX_prom_2 %>% dplyr::count(tx_name)
colnames(EBOX_transcripts_counts) <- c("transcript_id", "EBOX_n")

# Clean transcript IDs to remove version numbers (to match your data)
EBOX_transcripts_counts$transcript_id <- sub("\\..*", "", EBOX_transcripts_counts$transcript_id)

# Create comprehensive transcript annotation using your new data sources
all_transcript_ids <- unique(c(protein_coding$transcript_id, 
                               transcript_unique$transcript_id,
                               transcript_correlation_all$transcript_id))

transcript_annotation <- data.frame(
  transcript_id = all_transcript_ids,
  stringsAsFactors = FALSE
)

# Add EBOX counts
transcript_annotation <- left_join(transcript_annotation, EBOX_transcripts_counts, by = "transcript_id")
transcript_annotation$EBOX_n[is.na(transcript_annotation$EBOX_n)] <- 0

# Add group labels using your new classification system
transcript_annotation <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation", 
    transcript_id %in% protein_coding$transcript_id ~ "All",
    TRUE ~ "Other"
  )) %>%
  filter(group %in% c("All", "Correlation", "Unique"))

# Create EBOX distribution (same format as your original charts)
ebox_distribution <- transcript_annotation %>%
  dplyr::mutate(EBOX_n_cat = case_when(
    EBOX_n == 0 ~ "0",
    EBOX_n == 1 ~ "1",
    EBOX_n == 2 ~ "2",
    EBOX_n >= 3 ~ "≥3"
  )) %>%
  dplyr::count(group, EBOX_n_cat) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(total = sum(n), proportion = n / total) %>%
  dplyr::ungroup() %>%  # Add this line
  dplyr::select(group, EBOX_n_cat, n, proportion)

# Order the factor levels for plotting
ebox_distribution$EBOX_n_cat <- factor(ebox_distribution$EBOX_n_cat, 
                                       levels = c("0", "1", "2", "≥3"))

# Define custom colors (same as your other charts)
custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)

# Plot EBOX distribution (same formatting as your other charts)
ggplot(ebox_distribution, aes(x = EBOX_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  # labs(title = "Number of E-Box", y = "Percentage of Transcripts") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )




######### Statistical Significance:
library(dplyr)
library(tidyr)

# Step 1: Prepare categorized data
df_stats <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(EBOX_n_cat = case_when(
    EBOX_n == 0 ~ "0",
    EBOX_n == 1 ~ "1",
    EBOX_n == 2 ~ "2",
    EBOX_n >= 3 ~ "≥3"
  ))

# Step 2: Define groups and categories
groups <- c("Unique", "Correlation", "All")
categories <- c("0", "1", "2", "≥3")

# Step 3: Function to compute Fisher test between two groups for a category
fisher_compare <- function(group1, group2, cat) {
  data_sub <- df_stats %>%
    filter(group %in% c(group1, group2)) %>%
    mutate(in_category = EBOX_n_cat == cat)
  
  # Contingency table
  tab <- table(data_sub$group, data_sub$in_category)
  
  # Run Fisher's exact test with simulation
  test <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  
  data.frame(
    Comparison = paste(group1, "vs", group2),
    Category = cat,
    P_value = signif(test$p.value, 4),
    stringsAsFactors = FALSE
  )
}

# Step 4: Run all pairwise comparisons
results <- list()

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    for (cat in categories) {
      res <- fisher_compare(g1, g2, cat)
      results[[length(results) + 1]] <- res
    }
  }
}

# Step 5: Combine into a data frame and print
all_ebox_results <- do.call(rbind, results)

# Add multiple testing correction and enhanced output
all_ebox_results$P_adjusted_FDR <- p.adjust(all_ebox_results$P_value, method = "fdr")

# Add significance indicators
all_ebox_results$Significance_raw <- case_when(
  all_ebox_results$P_value < 0.0001 ~ "****",
  all_ebox_results$P_value < 0.001 ~ "***", 
  all_ebox_results$P_value < 0.01 ~ "**",
  all_ebox_results$P_value < 0.05 ~ "*",
  TRUE ~ "n.s."
)

all_ebox_results$Significance_adjusted <- case_when(
  all_ebox_results$P_adjusted_FDR < 0.0001 ~ "****",
  all_ebox_results$P_adjusted_FDR < 0.001 ~ "***",
  all_ebox_results$P_adjusted_FDR < 0.01 ~ "**", 
  all_ebox_results$P_adjusted_FDR < 0.05 ~ "*",
  TRUE ~ "n.s."
)

print(all_ebox_results)

write.csv(transcript_annotation, "output/transcript_ebox_counts.csv", row.names = FALSE)
