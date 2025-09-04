library(dplyr)
library(ggplot2)

PRJEB30337 <- read.csv("Data/PRJEB30337_TPM.csv")
GSE163646 <- read.csv("Data/GSE_163646_kallisto_transcript_TPM.csv")
Henja <- read.csv("Data/Henja_TPM.csv")
GSE283655 <- read.csv("Data/GSE283655_kallisto_transcript_TPM.csv")
# Trim version suffix from transcript_id in Henja and PRJEB30337
# Trim version suffix from transcript_id in Henja and PRJEB30337
Henja$transcript_id <- sub("\\.\\d+$", "", Henja$transcript_id)
PRJEB30337$transcript_id <- sub("\\.\\d+$", "", PRJEB30337$transcript_id)


# Sequential full merges
merged_siMITF <- merge(PRJEB30337, GSE163646, by = "transcript_id", all = TRUE)
merged_siMITF <- merge(merged_siMITF, Henja, by = "transcript_id", all = TRUE)
merged_siMITF <- merge(merged_siMITF, GSE283655, by = "transcript_id", all = TRUE)

# Define sample groups
PRJEB30337_con <- c("ERR3013912",	"ERR3013913")
PRJEB30337_siMITF  <- c("ERR3013908",	"ERR3013909")
GSE163646_CON    <- c("SRR13282348", "SRR13282350")
GSE163646_siMTIF   <- c("SRR13282352", "SRR13282353")
GSE283655_con <- c("SRR31631254",	"SRR31631258",	"SRR31631262")
GSE283655_siMITF  <- c("SRR31631252",	"SRR31631256", "SRR31631260")
Henja_con <- c("SSRR7346993", "SRR7346994")
Henja_siMITF  <- c("SRR7346991", "SRR7346992")

# Calculate mean TPM per group
mean_expr_df <- merged_siMITF %>%
  mutate(
    PRJEB30337_CON_mean     = rowMeans(select(., all_of(PRJEB30337_con)), na.rm = TRUE),
    PRJEB30337_siMITF_mean  = rowMeans(select(., all_of(PRJEB30337_siMITF)), na.rm = TRUE),
    GSE163646_CON_mean      = rowMeans(select(., all_of(GSE163646_CON)), na.rm = TRUE),
    GSE163646_siMITF_mean   = rowMeans(select(., all_of(GSE163646_siMITF)), na.rm = TRUE),
    GSE283655_CON_mean      = rowMeans(select(., all_of(GSE283655_con)), na.rm = TRUE),
    GSE283655_siMITF_mean   = rowMeans(select(., all_of(GSE283655_siMITF)), na.rm = TRUE),
    Henja_CON_mean          = rowMeans(select(., all_of(Henja_con)), na.rm = TRUE),
    Henja_siMITF_mean       = rowMeans(select(., all_of(Henja_siMITF)), na.rm = TRUE)
  ) %>%
  select(transcript_id,
         PRJEB30337_CON_mean, PRJEB30337_siMITF_mean,
         GSE163646_CON_mean,  GSE163646_siMITF_mean,
         GSE283655_CON_mean,  GSE283655_siMITF_mean,
         Henja_CON_mean,      Henja_siMITF_mean)

# Compute ratios
ratio_df <- mean_expr_df %>%
  mutate(
    PRJEB30337_ratio = PRJEB30337_siMITF_mean / PRJEB30337_CON_mean,
    GSE163646_ratio  = GSE163646_siMITF_mean  / GSE163646_CON_mean,
    GSE283655_ratio  = GSE283655_siMITF_mean  / GSE283655_CON_mean,
    Henja_ratio      = Henja_siMITF_mean      / Henja_CON_mean
  ) %>%
  select(transcript_id,
         PRJEB30337_ratio, GSE163646_ratio,
         GSE283655_ratio, Henja_ratio) %>%
  filter(if_all(-transcript_id, is.finite))

# Load groups
protein_coding <- read.csv("Data/protein_coding_RESUBMISSION.csv")
correlated     <- read.csv("Data/correlated_RESUBMISSION.csv")
discordant     <- read.csv("Data/Discordant_RESUBMISSION.csv")
library(dplyr)
correlated <- correlated %>%
  filter(!(transcript_id %in% discordant$transcript_id))

# Annotate group
ratio_df_annotated <- ratio_df %>%
  mutate(group = case_when(
    transcript_id %in% discordant$transcript_id     ~ "discordant",
    transcript_id %in% correlated$transcript_id     ~ "correlated",
    transcript_id %in% protein_coding$transcript_id ~ "protein_coding",
    TRUE                                             ~ "other"
  )) %>%
  filter(group %in% c("protein_coding", "correlated", "discordant"))

# Binning
ratio_df_binned <- ratio_df_annotated %>%
  mutate(
    PRJEB30337_bin = cut(PRJEB30337_ratio, breaks = c(-Inf, 0.5, 1, Inf), labels = c("≤0.5", "0.5–1", "≥1")),
    GSE163646_bin  = cut(GSE163646_ratio,  breaks = c(-Inf, 0.5, 1, Inf), labels = c("≤0.5", "0.5–1", "≥1")),
    GSE283655_bin  = cut(GSE283655_ratio,  breaks = c(-Inf, 0.5, 1, Inf), labels = c("≤0.5", "0.5–1", "≥1")),
    Henja_bin      = cut(Henja_ratio,      breaks = c(-Inf, 0.5, 1, Inf), labels = c("≤0.5", "0.5–1", "≥1"))
  )

# Plotting function
plot_ratio_bins <- function(df, dataset, bin_col, title_suffix) {
  df %>%
    count(group, .data[[bin_col]]) %>%
    ggplot(aes(x = group, y = n, fill = .data[[bin_col]])) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("siMITF / CON Ratio Bins", title_suffix),
      x = "Transcript Group",
      y = "Transcript Count",
      fill = "Ratio Bin"
    ) +
    theme_minimal()
}

# Generate all bar plots
plot_ratio_bins(ratio_df_binned, "PRJEB30337", "PRJEB30337_bin", "(PRJEB30337)")
plot_ratio_bins(ratio_df_binned, "GSE163646",  "GSE163646_bin",  "(GSE163646)")
plot_ratio_bins(ratio_df_binned, "GSE283655",  "GSE283655_bin",  "(GSE283655)")
plot_ratio_bins(ratio_df_binned, "Henja",      "Henja_bin",      "(Henja)")

# Percent plots
plot_ratio_percent <- function(df, bin_col, dataset) {
  df %>%
    count(group, .data[[bin_col]]) %>%
    group_by(group) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ggplot(aes(x = group, y = percent, fill = .data[[bin_col]])) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Percent of Transcripts by Ratio Bin", dataset),
      x = "Transcript Group",
      y = "Percent",
      fill = "Ratio Bin"
    ) +
    theme_minimal()
}

# Generate percent plots
plot_ratio_percent(ratio_df_binned, "PRJEB30337_bin", "PRJEB30337")
plot_ratio_percent(ratio_df_binned, "GSE163646_bin",  "GSE163646")
plot_ratio_percent(ratio_df_binned, "GSE283655_bin",  "GSE283655")
plot_ratio_percent(ratio_df_binned, "Henja_bin",      "Henja")







# Add NEW binning with 4 levels
ratio_df_binned4 <- ratio_df_annotated %>%
  mutate(
    PRJEB30337_bin4 = case_when(
      PRJEB30337_ratio <= 0.25 ~ "≤0.25",
      PRJEB30337_ratio <= 0.5  ~ "0.25–0.5",
      PRJEB30337_ratio <= 1    ~ "0.51–1",
      TRUE                     ~ "≥1"
    ),
    GSE163646_bin4 = case_when(
      GSE163646_ratio <= 0.25 ~ "≤0.25",
      GSE163646_ratio <= 0.5  ~ "0.25–0.5",
      GSE163646_ratio <= 1    ~ "0.51–1",
      TRUE                    ~ "≥1"
    ),
    GSE283655_bin4 = case_when(
      GSE283655_ratio <= 0.25 ~ "≤0.25",
      GSE283655_ratio <= 0.5  ~ "0.25–0.5",
      GSE283655_ratio <= 1    ~ "0.51–1",
      TRUE                    ~ "≥1"
    ),
    Henja_bin4 = case_when(
      Henja_ratio <= 0.25 ~ "≤0.25",
      Henja_ratio <= 0.5  ~ "0.25–0.5",
      Henja_ratio <= 1    ~ "0.51–1",
      TRUE                ~ "≥1"
    )
  )

# Function for stacked count bar plots (4 bins)
plot_ratio_bins_4 <- function(df, bin_col, title_suffix) {
  df %>%
    count(group, .data[[bin_col]]) %>%
    ggplot(aes(x = group, y = n, fill = .data[[bin_col]])) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("siMITF / CON 4-bin Ratio (", title_suffix, ")", sep = ""),
      x = "Transcript Group", y = "Transcript Count", fill = "Ratio Bin"
    ) +
    theme_minimal()
}

# Function for stacked percent bar plots (4 bins)
plot_ratio_percent_4 <- function(df, bin_col, title_suffix) {
  df %>%
    count(group, .data[[bin_col]]) %>%
    group_by(group) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ggplot(aes(x = group, y = percent, fill = .data[[bin_col]])) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Percent by siMITF / CON 4-bin Ratio (", title_suffix, ")", sep = ""),
      x = "Transcript Group", y = "Percent", fill = "Ratio Bin"
    ) +
    theme_minimal()
}

# 🔢 COUNT PLOTS
plot_ratio_bins_4(ratio_df_binned4, "PRJEB30337_bin4", "PRJEB30337")
plot_ratio_bins_4(ratio_df_binned4, "GSE163646_bin4",  "GSE163646")
plot_ratio_bins_4(ratio_df_binned4, "GSE283655_bin4",  "GSE283655")
plot_ratio_bins_4(ratio_df_binned4, "Henja_bin4",      "Henja")

# 📊 PERCENT PLOTS
plot_ratio_percent_4(ratio_df_binned4, "PRJEB30337_bin4", "PRJEB30337")
plot_ratio_percent_4(ratio_df_binned4, "GSE163646_bin4",  "GSE163646")
plot_ratio_percent_4(ratio_df_binned4, "GSE283655_bin4",  "GSE283655")
plot_ratio_percent_4(ratio_df_binned4, "Henja_bin4",      "Henja")




# Function to run pairwise Fisher's tests for a given bin column
pairwise_bin_tests <- function(df, bin_col, dataset_name) {
  bin_levels <- sort(unique(na.omit(df[[bin_col]])))
  
  group_pairs <- list(
    c("discordant", "correlated"),
    c("discordant", "protein_coding"),
    c("correlated", "protein_coding")
  )
  
  results <- list()
  
  for (bin in bin_levels) {
    for (pair in group_pairs) {
      g1 <- pair[1]
      g2 <- pair[2]
      
      sub_df <- df %>% filter(group %in% c(g1, g2))
      tab <- table(
        group = sub_df$group,
        in_bin = sub_df[[bin_col]] == bin
      )
      
      if (all(c(g1, g2) %in% rownames(tab))) {
        test <- fisher.test(tab)
        results[[length(results) + 1]] <- data.frame(
          dataset = dataset_name,
          binning = bin_col,
          bin = bin,
          group1 = g1,
          group2 = g2,
          p_value = test$p.value,
          odds_ratio = test$estimate,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  bind_rows(results)
}

# Run tests on 4-bin data
res4a <- pairwise_bin_tests(ratio_df_binned4, "PRJEB30337_bin4", "PRJEB30337")
res4b <- pairwise_bin_tests(ratio_df_binned4, "GSE163646_bin4",  "GSE163646")
res4c <- pairwise_bin_tests(ratio_df_binned4, "GSE283655_bin4",  "GSE283655")
res4d <- pairwise_bin_tests(ratio_df_binned4, "Henja_bin4",      "Henja")
results_4bin <- bind_rows(res4a, res4b, res4c, res4d)

# Run tests on 3-bin data
res3a <- pairwise_bin_tests(ratio_df_binned, "PRJEB30337_bin", "PRJEB30337")
res3b <- pairwise_bin_tests(ratio_df_binned, "GSE163646_bin",  "GSE163646")
res3c <- pairwise_bin_tests(ratio_df_binned, "GSE283655_bin",  "GSE283655")
res3d <- pairwise_bin_tests(ratio_df_binned, "Henja_bin",      "Henja")
results_3bin <- bind_rows(res3a, res3b, res3c, res3d)

# Combine all results and adjust p-values
all_results <- bind_rows(results_3bin, results_4bin) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"),
         significance = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE          ~ ""
         ))

# View top results
head(all_results)

# Optional: save to CSV
# write.csv(all_results, "pairwise_bin_fisher_results.csv", row.names = FALSE)




library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Compute log2 ratios
log2_df <- ratio_df_annotated %>%
  mutate(
    PRJEB30337_log2 = log2(PRJEB30337_ratio),
    GSE163646_log2  = log2(GSE163646_ratio),
    GSE283655_log2  = log2(GSE283655_ratio),
    Henja_log2      = log2(Henja_ratio)
  ) %>%
  select(transcript_id, group,
         PRJEB30337_log2, GSE163646_log2, GSE283655_log2, Henja_log2)

# Step 2: Pivot to long format
log2_long <- log2_df %>%
  pivot_longer(
    cols = ends_with("_log2"),
    names_to = "dataset",
    values_to = "log2_ratio"
  ) %>%
  filter(group %in% c("protein_coding", "correlated", "discordant")) %>%
  filter(is.finite(log2_ratio))

# Step 3: Factor ordering and custom colors
log2_long$group <- factor(log2_long$group,
                          levels = c("protein_coding", "correlated", "discordant"))

group_colors <- c(
  "protein_coding" = "lightpink",
  "correlated" = "hotpink",
  "discordant" = "red"
)

# Step 4: Summary stats for plotting
summary_stats <- log2_long %>%
  group_by(group, dataset) %>%
  summarise(
    mean_log2 = mean(log2_ratio, na.rm = TRUE),
    sem_log2  = sd(log2_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Clean dataset names
summary_stats$dataset <- gsub("_log2", "", summary_stats$dataset)
log2_long$dataset <- gsub("_log2", "", log2_long$dataset)

# Step 5: Compute pairwise p-values for significance
pairwise_pvals <- log2_long %>%
  group_by(dataset) %>%
  do({
    comparisons <- list(c("protein_coding", "correlated"),
                        c("protein_coding", "discordant"),
                        c("correlated", "discordant"))
    
    res <- lapply(comparisons, function(pair) {
      test <- t.test(log2_ratio ~ group,
                     data = filter(., group %in% pair))
      data.frame(
        dataset = unique(.$dataset),
        group1 = pair[1],
        group2 = pair[2],
        p = test$p.value
      )
    })
    bind_rows(res)
  }) %>%
  mutate(star = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  ))

# Step 6: Plot with bar + SEM + significance
ggplot(summary_stats, aes(x = group, y = mean_log2, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_log2 - sem_log2, ymax = mean_log2 + sem_log2),
                width = 0.2) +
  facet_wrap(~dataset) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Average log2(siMITF / CON) per Transcript Group",
    x = "Transcript Group",
    y = "log2(siMITF / CON)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold")) +
  geom_text(data = pairwise_pvals,
            aes(
              x = (as.numeric(factor(group1, levels = levels(summary_stats$group))) +
                     as.numeric(factor(group2, levels = levels(summary_stats$group))) ) / 2,
              y = max(summary_stats$mean_log2) + 0.5,
              label = star
            ),
            inherit.aes = FALSE,
            size = 5)
