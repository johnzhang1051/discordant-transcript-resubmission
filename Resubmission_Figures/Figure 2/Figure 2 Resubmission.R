library(dplyr)
library(ggplot2)
library(tidyr)
library(textshape)

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
  dplyr::select(group, MITF_peak_n_cat, n, proportion)

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
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  #labs(title = "Number of MITF ChIP Peaks",y = "Percentage of Transcripts") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


## Average Knockdown Distribution (same formatting as original)
# Clean avg_knockdown: remove NaN, zeros, and values > 2
knockdown_data <- transcript_annotation %>%
  filter(!is.na(avg_knockdown) & !is.nan(avg_knockdown)) %>%
  filter(avg_knockdown > 0 & avg_knockdown <= 2) %>%
  dplyr::select(transcript_id, group, avg_knockdown)

# Boxplot (same formatting as original)
ggplot(knockdown_data, aes(x = group, y = avg_knockdown, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  #labs(title = "Transcript Expression (siMITF/siCON)",y = "Relative Expression") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


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
  dplyr::select(group, knockdown_bin, n, proportion)

# Set the order of knockdown_bin
knockdown_binned$knockdown_bin <- factor(knockdown_binned$knockdown_bin,
                                         levels = c("≤0.25", "0.26–0.5", "0.51–1", "≥1"))

# Plot binned knockdown (same formatting as original)
ggplot(knockdown_binned, aes(x = knockdown_bin, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  #labs(title = "Transcript Expression (siMITF/siCON) (binned)", y = "Percentage of Transcripts") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.position = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())



## Proportion of Transcripts with Unique Promoters
library(GenomicRanges)

# Extract transcript coordinates from Kenny dataset for promoter analysis
transcript_coords <- Kenny %>%
  dplyr::select(transcript_id, seqnames, start, end, strand, geneId) %>%
  distinct() %>%
  mutate(
    # Define promoter as 1kb upstream of TSS
    promoter_start = ifelse(strand == 1, start - 1000, end),
    promoter_end = ifelse(strand == 1, start, end + 1000),
    promoter_start = pmax(1, promoter_start)
  )

transcript_coords_clean <- transcript_coords %>%
  dplyr::select(-start, -end)  # Remove conflicting column names

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
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

### ### ### ### ### ### STATISTICAL TESTS FOR P-VALUES ### ### ### ### ### ### 

library(dplyr)
library(tidyr)

# 1. MITF PEAK ANALYSIS - Fisher's exact tests
cat("=== MITF PEAK ANALYSIS ===\n")

# Prepare MITF data for testing
df_stats_mitf <- transcript_annotation %>%
  mutate(MITF_peak_n_cat = case_when(
    MITF_peak_n == 0 ~ "0",
    MITF_peak_n == 1 ~ "1",
    MITF_peak_n >= 2 ~ "≥2"
  ))

# Define groups and categories
groups <- c("Unique", "Correlation", "All")
mitf_categories <- c("0", "1", "≥2")

# Function for Fisher test on MITF peaks
fisher_compare_mitf <- function(group1, group2, cat) {
  data_sub <- df_stats_mitf %>%
    filter(group %in% c(group1, group2)) %>%
    mutate(in_category = MITF_peak_n_cat == cat)
  
  tab <- table(data_sub$group, data_sub$in_category)
  test <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  
  data.frame(
    Comparison = paste(group1, "vs", group2),
    Category = cat,
    P_value = signif(test$p.value, 4),
    Odds_Ratio = as.numeric(test$estimate),
    CI_lower = test$conf.int[1],
    CI_upper = test$conf.int[2],
    stringsAsFactors = FALSE
  )
}

# Run MITF comparisons
mitf_results <- list()
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    for (cat in mitf_categories) {
      res <- fisher_compare_mitf(g1, g2, cat)
      mitf_results[[length(mitf_results) + 1]] <- res
    }
  }
}

mitf_all_results <- do.call(rbind, mitf_results)
mitf_all_results$P_adjusted_FDR <- p.adjust(mitf_all_results$P_value, method = "fdr")

# Add significance indicators
mitf_all_results$Significance_raw <- case_when(
  mitf_all_results$P_value < 0.0001 ~ "****",
  mitf_all_results$P_value < 0.001 ~ "***",
  mitf_all_results$P_value < 0.01 ~ "**",
  mitf_all_results$P_value < 0.05 ~ "*",
  TRUE ~ "n.s."
)

mitf_all_results$Significance_adjusted <- case_when(
  mitf_all_results$P_adjusted_FDR < 0.0001 ~ "****",
  mitf_all_results$P_adjusted_FDR < 0.001 ~ "***",
  mitf_all_results$P_adjusted_FDR < 0.01 ~ "**",
  mitf_all_results$P_adjusted_FDR < 0.05 ~ "*",
  TRUE ~ "n.s."
)

print("MITF Peak Analysis Results:")
print(mitf_all_results)


# 2. AVERAGE KNOCKDOWN ANALYSIS - Wilcoxon tests for continuous data
cat("\n=== AVERAGE KNOCKDOWN ANALYSIS (Continuous) ===\n")

# Use the same knockdown_data you already created
knockdown_stats_results <- list()

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    
    group1_data <- knockdown_data %>% filter(group == g1) %>% pull(avg_knockdown)
    group2_data <- knockdown_data %>% filter(group == g2) %>% pull(avg_knockdown)
    
    # Wilcoxon test
    wilcox_test <- wilcox.test(group1_data, group2_data)
    
    # Effect size (r = Z / sqrt(N))
    n1 <- length(group1_data)
    n2 <- length(group2_data)
    n_total <- n1 + n2
    r_effect_size <- abs(qnorm(wilcox_test$p.value/2)) / sqrt(n_total)
    
    effect_interpretation <- case_when(
      r_effect_size < 0.1 ~ "negligible",
      r_effect_size < 0.3 ~ "small",
      r_effect_size < 0.5 ~ "medium",
      TRUE ~ "large"
    )
    
    knockdown_stats_results[[length(knockdown_stats_results) + 1]] <- data.frame(
      Comparison = paste(g1, "vs", g2),
      P_value = wilcox_test$p.value,
      Effect_size_r = r_effect_size,
      Effect_interpretation = effect_interpretation,
      Group1_median = median(group1_data, na.rm = TRUE),
      Group2_median = median(group2_data, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

knockdown_continuous_results <- do.call(rbind, knockdown_stats_results)
knockdown_continuous_results$P_adjusted_FDR <- p.adjust(knockdown_continuous_results$P_value, method = "fdr")

knockdown_continuous_results$Significance_raw <- case_when(
  knockdown_continuous_results$P_value < 0.0001 ~ "****",
  knockdown_continuous_results$P_value < 0.001 ~ "***",
  knockdown_continuous_results$P_value < 0.01 ~ "**",
  knockdown_continuous_results$P_value < 0.05 ~ "*",
  TRUE ~ "n.s."
)

knockdown_continuous_results$Significance_adjusted <- case_when(
  knockdown_continuous_results$P_adjusted_FDR < 0.0001 ~ "****",
  knockdown_continuous_results$P_adjusted_FDR < 0.001 ~ "***",
  knockdown_continuous_results$P_adjusted_FDR < 0.01 ~ "**",
  knockdown_continuous_results$P_adjusted_FDR < 0.05 ~ "*",
  TRUE ~ "n.s."
)

print("Average Knockdown Analysis Results:")
print(knockdown_continuous_results)


# 3. BINNED KNOCKDOWN ANALYSIS - Fisher's exact tests
cat("\n=== BINNED KNOCKDOWN ANALYSIS ===\n")

# Prepare binned knockdown data for testing
df_stats_knockdown_binned <- transcript_annotation %>%
  filter(!is.na(avg_knockdown) & !is.nan(avg_knockdown), avg_knockdown > 0) %>%
  mutate(knockdown_bin = case_when(
    avg_knockdown <= 0.25 ~ "≤0.25",
    avg_knockdown <= 0.5 ~ "0.26–0.5",
    avg_knockdown <= 1 ~ "0.51–1",
    avg_knockdown > 1 ~ "≥1"
  ))

knockdown_bins <- c("≤0.25", "0.26–0.5", "0.51–1", "≥1")

# Function for Fisher test on binned knockdown
fisher_compare_knockdown_binned <- function(group1, group2, bin) {
  data_sub <- df_stats_knockdown_binned %>%
    filter(group %in% c(group1, group2)) %>%
    mutate(in_bin = knockdown_bin == bin)
  
  tab <- table(data_sub$group, data_sub$in_bin)
  test <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  
  data.frame(
    Comparison = paste(group1, "vs", group2),
    Knockdown_Bin = bin,
    P_value = signif(test$p.value, 4),
    Odds_Ratio = as.numeric(test$estimate),
    CI_lower = test$conf.int[1],
    CI_upper = test$conf.int[2],
    stringsAsFactors = FALSE
  )
}

# Run binned knockdown comparisons
knockdown_binned_results <- list()
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    for (bin in knockdown_bins) {
      res <- fisher_compare_knockdown_binned(g1, g2, bin)
      knockdown_binned_results[[length(knockdown_binned_results) + 1]] <- res
    }
  }
}

knockdown_binned_all_results <- do.call(rbind, knockdown_binned_results)
knockdown_binned_all_results$P_adjusted_FDR <- p.adjust(knockdown_binned_all_results$P_value, method = "fdr")

knockdown_binned_all_results$Significance_raw <- case_when(
  knockdown_binned_all_results$P_value < 0.0001 ~ "****",
  knockdown_binned_all_results$P_value < 0.001 ~ "***",
  knockdown_binned_all_results$P_value < 0.01 ~ "**",
  knockdown_binned_all_results$P_value < 0.05 ~ "*",
  TRUE ~ "n.s."
)

knockdown_binned_all_results$Significance_adjusted <- case_when(
  knockdown_binned_all_results$P_adjusted_FDR < 0.0001 ~ "****",
  knockdown_binned_all_results$P_adjusted_FDR < 0.001 ~ "***",
  knockdown_binned_all_results$P_adjusted_FDR < 0.01 ~ "**",
  knockdown_binned_all_results$P_adjusted_FDR < 0.05 ~ "*",
  TRUE ~ "n.s."
)

print("Binned Knockdown Analysis Results:")
print(knockdown_binned_all_results)


# 4. UNIQUE PROMOTER ANALYSIS - Fisher's exact tests
cat("\n=== UNIQUE PROMOTER ANALYSIS ===\n")

# Function for unique promoter Fisher test
fisher_compare_unique_promoter <- function(group1, group2, transcript_list1, transcript_list2, unique_promoter_transcripts) {
  
  # Count overlaps for each group
  group1_unique <- sum(transcript_list1$transcript_id %in% unique_promoter_transcripts)
  group1_total <- nrow(transcript_list1)
  group1_non_unique <- group1_total - group1_unique
  
  group2_unique <- sum(transcript_list2$transcript_id %in% unique_promoter_transcripts)
  group2_total <- nrow(transcript_list2)
  group2_non_unique <- group2_total - group2_unique
  
  # Create contingency table
  tab <- matrix(c(group1_unique, group1_non_unique, group2_unique, group2_non_unique), 
                nrow = 2, 
                dimnames = list(c(group1, group2), c("Unique", "Non-unique")))
  
  test <- fisher.test(tab)
  
  data.frame(
    Comparison = paste(group1, "vs", group2),
    P_value = signif(test$p.value, 4),
    Odds_Ratio = as.numeric(test$estimate),
    CI_lower = test$conf.int[1],
    CI_upper = test$conf.int[2],
    Group1_proportion = group1_unique / group1_total,
    Group2_proportion = group2_unique / group2_total,
    stringsAsFactors = FALSE
  )
}

# Run unique promoter comparisons
unique_promoter_results <- list()

# All vs Correlation
res1 <- fisher_compare_unique_promoter("All", "Correlation", protein_coding, transcript_correlation_all, unique_promoter_transcripts)
unique_promoter_results[[1]] <- res1

# All vs Unique  
res2 <- fisher_compare_unique_promoter("All", "Unique", protein_coding, transcript_unique, unique_promoter_transcripts)
unique_promoter_results[[2]] <- res2

# Correlation vs Unique
res3 <- fisher_compare_unique_promoter("Correlation", "Unique", transcript_correlation_all, transcript_unique, unique_promoter_transcripts)
unique_promoter_results[[3]] <- res3

unique_promoter_all_results <- do.call(rbind, unique_promoter_results)
unique_promoter_all_results$P_adjusted_FDR <- p.adjust(unique_promoter_all_results$P_value, method = "fdr")

unique_promoter_all_results$Significance_raw <- case_when(
  unique_promoter_all_results$P_value < 0.0001 ~ "****",
  unique_promoter_all_results$P_value < 0.001 ~ "***",
  unique_promoter_all_results$P_value < 0.01 ~ "**",
  unique_promoter_all_results$P_value < 0.05 ~ "*",
  TRUE ~ "n.s."
)

unique_promoter_all_results$Significance_adjusted <- case_when(
  unique_promoter_all_results$P_adjusted_FDR < 0.0001 ~ "****",
  unique_promoter_all_results$P_adjusted_FDR < 0.001 ~ "***",
  unique_promoter_all_results$P_adjusted_FDR < 0.01 ~ "**",
  unique_promoter_all_results$P_adjusted_FDR < 0.05 ~ "*",
  TRUE ~ "n.s."
)

print("Unique Promoter Analysis Results:")
print(unique_promoter_all_results)


# 5. SUMMARY OF ALL SIGNIFICANT FINDINGS

# Collect all significant results
sig_mitf <- mitf_all_results %>% filter(P_adjusted_FDR < 0.05) %>% 
  mutate(Analysis = "MITF Peaks", Variable = Category)

sig_knockdown_cont <- knockdown_continuous_results %>% filter(P_adjusted_FDR < 0.05) %>%
  mutate(Analysis = "Knockdown (Continuous)", Variable = "avg_knockdown")

sig_knockdown_binned <- knockdown_binned_all_results %>% filter(P_adjusted_FDR < 0.05) %>%
  mutate(Analysis = "Knockdown (Binned)", Variable = Knockdown_Bin)

sig_unique_promoter <- unique_promoter_all_results %>% filter(P_adjusted_FDR < 0.05) %>%
  mutate(Analysis = "Unique Promoter", Variable = "unique_promoter")

# Print summaries
if(nrow(sig_mitf) > 0) {
  cat("Significant MITF Peak differences:\n")
  print(sig_mitf %>% dplyr::select(Analysis, Comparison, Variable, P_value, Significance_raw, P_adjusted_FDR, Significance_adjusted))
  cat("\n")
}

if(nrow(sig_knockdown_cont) > 0) {
  cat("Significant Knockdown (Continuous) differences:\n")
  print(sig_knockdown_cont %>% dplyr::select(Analysis, Comparison, P_value, Significance_raw, P_adjusted_FDR, Significance_adjusted, Effect_interpretation))
  cat("\n")
}

if(nrow(sig_knockdown_binned) > 0) {
  cat("Significant Knockdown (Binned) differences:\n")
  print(sig_knockdown_binned %>% dplyr::select(Analysis, Comparison, Variable, P_value, Significance_raw, P_adjusted_FDR, Significance_adjusted))
  cat("\n")
}

if(nrow(sig_unique_promoter) > 0) {
  cat("Significant Unique Promoter differences:\n")
  print(sig_unique_promoter %>% dplyr::select(Analysis, Comparison, P_value, Significance_raw, P_adjusted_FDR, Significance_adjusted))
  cat("\n")
}
