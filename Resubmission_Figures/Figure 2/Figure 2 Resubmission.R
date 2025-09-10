library(dplyr)
library(ggplot2)
library(tidyr)

# Load data using the same sources as your first two codes
protein_coding <- read.csv("data/protein_coding_RESUBMISSION.csv")
transcript_unique <- read.csv("data/discordant_RESUBMISSION.csv")  # discordant transcripts
transcript_correlation_all <- read.csv("data/correlated_RESUBMISSION.csv")  # correlated transcripts

# Load ChIP-seq data (from first code)
Kenny <- read.csv("data/Kenny.csv")
Laurette <- read.csv("data/Laurette.csv")
Louph <- read.csv("data/Louphrasitthiphol.csv")

# Standardize transcript IDs in ChIP-seq data
standardize_transcripts <- function(df) {
  df %>%
    dplyr::rename(transcript_id = transcriptId) %>%
    dplyr::mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))
}

Kenny <- standardize_transcripts(Kenny)
Laurette <- standardize_transcripts(Laurette)
Louph <- standardize_transcripts(Louph)

# Filter for promoter peaks (<=1kb)
filter_promoters <- function(df) {
  df %>%
    dplyr::filter(annotation == "Promoter (<=1kb)") %>%
    dplyr::select(transcript_id, annotation)
}

Kenny_promoters <- filter_promoters(Kenny)
Laurette_promoters <- filter_promoters(Laurette)
Louph_promoters <- filter_promoters(Louph)

# Load knockdown data (from second code)
PRJEB30337 <- read.csv("data/PRJEB30337_TPM.csv")
GSE163646 <- read.csv("data/GSE_163646_kallisto_transcript_TPM.csv")
Henja <- read.csv("data/Henja_TPM.csv")
GSE283655 <- read.csv("data/GSE283655_kallisto_transcript_TPM.csv")

# Clean transcript IDs
Henja$transcript_id <- sub("\\.\\d+$", "", Henja$transcript_id)
PRJEB30337$transcript_id <- sub("\\.\\d+$", "", PRJEB30337$transcript_id)

# Merge all knockdown datasets
merged_siMITF <- Reduce(function(x, y) merge(x, y, by = "transcript_id", all = TRUE),
                        list(PRJEB30337, GSE163646, Henja, GSE283655))

# Define sample groups for knockdown analysis
PRJEB30337_con <- c("ERR3013912", "ERR3013913")
PRJEB30337_siMITF <- c("ERR3013908", "ERR3013909")
GSE163646_CON <- c("SRR13282348", "SRR13282350")
GSE163646_siMITF <- c("SRR13282352", "SRR13282353")
GSE283655_con <- c("SRR31631254", "SRR31631258", "SRR31631262")
GSE283655_siMITF <- c("SRR31631252", "SRR31631256", "SRR31631260")
Henja_con <- c("SRR7346993", "SRR7346994")
Henja_siMITF <- c("SRR7346991", "SRR7346992")

# Calculate mean TPM per group
mean_expr_df <- merged_siMITF %>%
  mutate(
    PRJEB30337_CON_mean = rowMeans(dplyr::select(., all_of(PRJEB30337_con)), na.rm = TRUE),
    PRJEB30337_siMITF_mean = rowMeans(dplyr::select(., all_of(PRJEB30337_siMITF)), na.rm = TRUE),
    GSE163646_CON_mean = rowMeans(dplyr::select(., all_of(GSE163646_CON)), na.rm = TRUE),
    GSE163646_siMITF_mean = rowMeans(dplyr::select(., all_of(GSE163646_siMITF)), na.rm = TRUE),
    GSE283655_CON_mean = rowMeans(dplyr::select(., all_of(GSE283655_con)), na.rm = TRUE),
    GSE283655_siMITF_mean = rowMeans(dplyr::select(., all_of(GSE283655_siMITF)), na.rm = TRUE),
    Henja_CON_mean = rowMeans(dplyr::select(., all_of(Henja_con)), na.rm = TRUE),
    Henja_siMITF_mean = rowMeans(dplyr::select(., all_of(Henja_siMITF)), na.rm = TRUE)
  ) %>%
  dplyr::select(transcript_id, PRJEB30337_CON_mean, PRJEB30337_siMITF_mean,
                GSE163646_CON_mean, GSE163646_siMITF_mean,
                GSE283655_CON_mean, GSE283655_siMITF_mean,
                Henja_CON_mean, Henja_siMITF_mean)

# Compute knockdown ratios
ratio_df <- mean_expr_df %>%
  mutate(
    PRJEB30337_ratio = PRJEB30337_siMITF_mean / PRJEB30337_CON_mean,
    GSE163646_ratio = GSE163646_siMITF_mean / GSE163646_CON_mean,
    GSE283655_ratio = GSE283655_siMITF_mean / GSE283655_CON_mean,
    Henja_ratio = Henja_siMITF_mean / Henja_CON_mean
  ) %>%
  dplyr::select(transcript_id, PRJEB30337_ratio, GSE163646_ratio,
                GSE283655_ratio, Henja_ratio) %>%
  filter(if_all(-transcript_id, is.finite))

# Calculate average knockdown across datasets
avg_knockdown_df <- ratio_df %>%
  mutate(avg_knockdown = rowMeans(dplyr::select(., -transcript_id), na.rm = TRUE)) %>%
  dplyr::select(transcript_id, avg_knockdown)

# Create comprehensive annotation dataset
transcript_annotation <- expand.grid(
  transcript_id = unique(c(protein_coding$transcript_id, 
                           transcript_unique$transcript_id,
                           transcript_correlation_all$transcript_id)),
  stringsAsFactors = FALSE
)

# Add ChIP-seq peak counts
transcript_annotation <- transcript_annotation %>%
  mutate(
    kenny_peak = as.integer(transcript_id %in% Kenny_promoters$transcript_id),
    laurette_peak = as.integer(transcript_id %in% Laurette_promoters$transcript_id),
    louph_peak = as.integer(transcript_id %in% Louph_promoters$transcript_id),
    MITF_peak_n = kenny_peak + laurette_peak + louph_peak
  )

# Add average knockdown data
transcript_annotation <- left_join(transcript_annotation, avg_knockdown_df, by = "transcript_id")

# Add group labels using the same logic as your first two codes
transcript_annotation <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation", 
    transcript_id %in% protein_coding$transcript_id ~ "All",
    TRUE ~ "Other"
  )) %>%
  filter(group %in% c("All", "Correlation", "Unique"))  # Keep only the main groups

## MITF Peak Distribution (same formatting as original)
mitf_distribution <- transcript_annotation %>%
  mutate(MITF_peak_n_cat = case_when(
    MITF_peak_n == 0 ~ "0",
    MITF_peak_n == 1 ~ "1",
    MITF_peak_n >= 2 ~ "≥2"
  )) %>%
  group_by(group, MITF_peak_n_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  select(group, MITF_peak_n_cat, n, proportion)

# Reorder factor levels for plotting
mitf_distribution$MITF_peak_n_cat <- factor(mitf_distribution$MITF_peak_n_cat,
                                            levels = c("0", "1", "≥2"))

# Define custom colors (same as original)
custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)

# Plot MITF_peak_n category proportions (same formatting as original)
ggplot(mitf_distribution, aes(x = MITF_peak_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Number of MITF ChIP Peaks",
       y = "Percentage of Transcripts") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(r = 15))
  )

## Average Knockdown Distribution (same formatting as original)
# Clean avg_knockdown: remove NaN, zeros, and values > 2
knockdown_data <- transcript_annotation %>%
  filter(!is.na(avg_knockdown) & !is.nan(avg_knockdown)) %>%
  filter(avg_knockdown > 0 & avg_knockdown <= 2) %>%
  select(transcript_id, group, avg_knockdown)

# Boxplot (same formatting as original)
ggplot(knockdown_data, aes(x = group, y = avg_knockdown, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Transcript Expression (siMITF/siCON)",
    y = "Relative Expression"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(r = 15))
  )

## Average Knockdown Binned Distribution (same formatting as original)
# Clean and bin avg_knockdown
knockdown_binned <- transcript_annotation %>%
  filter(!is.na(avg_knockdown) & !is.nan(avg_knockdown), avg_knockdown > 0) %>%
  mutate(knockdown_bin = case_when(
    avg_knockdown <= 0.25 ~ "≤0.25",
    avg_knockdown <= 0.5 ~ "0.26–0.5",
    avg_knockdown <= 1 ~ "0.51–1",
    avg_knockdown > 1 ~ "≥1"
  )) %>%
  group_by(group, knockdown_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  select(group, knockdown_bin, n, proportion)

# Set the order of knockdown_bin
knockdown_binned$knockdown_bin <- factor(knockdown_binned$knockdown_bin,
                                         levels = c("≤0.25", "0.26–0.5", "0.51–1", "≥1"))

# Plot binned knockdown (same formatting as original)
ggplot(knockdown_binned, aes(x = knockdown_bin, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Transcript Expression (siMITF/siCON) (binned)", y = "Percentage of Transcripts") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.title = element_blank())


## Proportion of Transcripts with Unique Promoters
library(GenomicRanges)

# Extract transcript coordinates from Kenny dataset for promoter analysis
transcript_coords <- Kenny %>%
  select(transcript_id, seqnames, start, end, strand, geneId) %>%
  distinct() %>%
  mutate(
    # Define promoter as 1kb upstream of TSS
    promoter_start = ifelse(strand == 1, start - 1000, end),
    promoter_end = ifelse(strand == 1, start, end + 1000),
    promoter_start = pmax(1, promoter_start)
  )

transcript_coords_clean <- transcript_coords %>%
  select(-start, -end)  # Remove conflicting column names

promoter_gr <- makeGRangesFromDataFrame(
  transcript_coords_clean,
  seqnames.field = "seqnames",
  start.field = "promoter_start", 
  end.field = "promoter_end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

# Find overlaps
overlaps <- findOverlaps(promoter_gr, promoter_gr)
overlap_summary <- as.data.frame(overlaps) %>%
  filter(queryHits != subjectHits) %>%
  group_by(queryHits) %>%
  summarise(n_overlaps = n(), .groups = "drop")

# Transcripts with unique (non-overlapping) promoters
unique_promoter_indices <- setdiff(1:length(promoter_gr), overlap_summary$queryHits)
unique_promoter_transcripts <- transcript_coords$transcript_id[unique_promoter_indices]

# Function to calculate unique promoter proportions
compute_unique_promoter_proportion <- function(transcript_group, label) {
  overlap_count <- sum(transcript_group$transcript_id %in% unique_promoter_transcripts)
  total <- nrow(transcript_group)
  
  data.frame(
    group = label,
    with_unique_promoter = overlap_count,
    total = total,
    proportion = overlap_count / total
  )
}

# Calculate proportions using updated group names
unique_promoter_stats <- bind_rows(
  compute_unique_promoter_proportion(protein_coding, "All"),
  compute_unique_promoter_proportion(transcript_correlation_all, "Correlation"),
  compute_unique_promoter_proportion(transcript_unique, "Unique")
)

unique_promoter_stats$group <- factor(unique_promoter_stats$group, 
                                      levels = c("All", "Correlation", "Unique"))

# Plot (same style as original with custom colors)
ggplot(unique_promoter_stats, aes(x = group, y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Proportion of Transcripts with Unique Promoters",
    x = "Transcript Group",
    y = "Proportion of Transcripts"
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))



