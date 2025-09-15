transcript_annotation <- readRDS("hg38_annotated_transcripts_EXPANDING.rds")

transcript_unique <- read.csv("Data/final_list_transcript_id.txt")
transcript_correlation_all <-read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL.csv")

library(dplyr)
library(ggplot2)
library(tidyr)

# Add a label to each transcript
transcript_annotation <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  ))




## proportion tables 
mitf_distribution <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
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

###
ebox_distribution <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(EBOX_n_cat = case_when(
    EBOX_n == 0 ~ "0",
    EBOX_n == 1 ~ "1",
    EBOX_n == 2 ~ "2",
    EBOX_n >= 3 ~ ">3"
  )) %>%
  group_by(group, EBOX_n_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  dplyr::select(group, EBOX_n_cat, n, proportion)

# Order the factor levels for plotting
ebox_distribution$EBOX_n_cat <- factor(ebox_distribution$EBOX_n_cat, 
                                       levels = c("0", "1", "2", ">3"))



library(ggplot2)
library(ggpubr)

custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)

# Plot MITF_peak_n category proportions
ggplot(mitf_distribution, aes(x = MITF_peak_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Proportion of Transcripts by MITF_peak_n Category",
       y = "Proportion") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(r = 15))  # larger + moved left
  )


ggplot(ebox_distribution, aes(x = EBOX_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Proportion of Transcripts by EBOX_n Category",
       y = "Proportion") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(r = 15))  # larger + moved left
  )






# OPTIONAL: Statistics (Chi-square on full contingency table)
mitf_chisq <- mitf_distribution %>%
  tidyr::pivot_wider(names_from = MITF_peak_n_cat, values_from = n) %>%
  column_to_rownames("group") %>%
  chisq.test()
print(mitf_chisq)

# Plot EBOX_n category proportions
ggplot(ebox_distribution, aes(x = EBOX_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Proportion of Transcripts by EBOX_n Category",
       x = "EBOX_n Category", y = "Proportion") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.title = element_blank())

### new ebox binning
ebox_distribution <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(EBOX_n_cat = case_when(
    EBOX_n == 0 ~ "0",
    EBOX_n == 1 ~ "1",
    EBOX_n == 2 ~ "2",
    EBOX_n == 3 ~ "3",
    EBOX_n == 4 ~ "4",
    EBOX_n >= 5 ~ "≥5"
  )) %>%
  group_by(group, EBOX_n_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  dplyr::select(group, EBOX_n_cat, n, proportion)
# Ensure EBOX_n_cat is ordered
ebox_distribution$EBOX_n_cat <- factor(ebox_distribution$EBOX_n_cat,
                                       levels = c("0", "1", "2", "3", "4", "≥5"))

custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)


ggplot(ebox_distribution, aes(x = EBOX_n_cat, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Proportion of Transcripts by EBOX_n Category",
       x = "EBOX_n Category", y = "Proportion") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.title = element_blank())







### AVG Knockdown Violin plot
# Clean avg_knockdown: remove NaN, zeros, and values > 2
knockdown_data <- transcript_annotation %>%
  filter(!is.nan(avg_knockdown)) %>%
  filter(avg_knockdown > 0 & avg_knockdown <= 2) %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  dplyr::select(transcript_id, group, avg_knockdown)

custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)
ggplot(knockdown_data, aes(x = group, y = avg_knockdown, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Filtered avg_knockdown per Transcript Group (0 < x ≤ 2)",
    y = "Average Knockdown"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_text(size = 24, margin = margin(r = 15))  # y-axis title styling
  )



##AVG Knockdown
###Clean and Bin
# Clean and bin avg_knockdown
knockdown_binned <- transcript_annotation %>%
  filter(!is.nan(avg_knockdown), avg_knockdown > 0) %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(knockdown_bin = case_when(
    avg_knockdown <= 0.25 ~ "<=0.25",
    avg_knockdown <= 0.5 ~ "0.26–0.5",
    avg_knockdown <= 0.75 ~ "0.51–0.75",
    avg_knockdown <= 1 ~ "0.76–1",
    avg_knockdown > 1 ~ ">1"
  )) %>%
  group_by(group, knockdown_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  dplyr::select(group, knockdown_bin, n, proportion)

# Set the order of knockdown_bin
knockdown_binned$knockdown_bin <- factor(knockdown_binned$knockdown_bin,
                                         levels = c("<=0.25", "0.26–0.5", "0.51–0.75", "0.76–1", ">1"))

# Define custom fill colors for groups
custom_colors <- c(
  "All" = "#FADADD",          # light pink
  "Correlation" = "#F08080",  # medium pink
  "Unique" = "#FF0000"        # red
)

ggplot(knockdown_binned, aes(x = knockdown_bin, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Proportion of avg_knockdown Bins per Transcript Group",
       x = "avg_knockdown Bin", y = "Proportion") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.title = element_blank())

### all_MITF_peak_ statistics
library(dplyr)
library(tidyr)

# Step 1: Prepare categorized data
df_stats <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(MITF_peak_n_cat = case_when(
    MITF_peak_n == 0 ~ "0",
    MITF_peak_n == 1 ~ "1",
    MITF_peak_n >= 2 ~ "≥2"
  ))

# Step 2: Define groups and categories
groups <- c("Unique", "Correlation", "All")
categories <- c("0", "1", "≥2")

# Step 3: Function to compute Fisher test between two groups for a category
fisher_compare <- function(group1, group2, cat) {
  data_sub <- df_stats %>%
    filter(group %in% c(group1, group2)) %>%
    mutate(in_category = MITF_peak_n_cat == cat)
  
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
all_MITF_results <- do.call(rbind, results)


###. EBOX comparisons
library(dplyr)
library(tidyr)

# Step 1: Prepare EBOX-categorized data
df_stats_ebox <- transcript_annotation %>%
  mutate(group = case_when(
    transcript_id %in% transcript_unique$transcript_id ~ "Unique",
    transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
    TRUE ~ "All"
  )) %>%
  mutate(EBOX_n_cat = case_when(
    EBOX_n == 0 ~ "0",
    EBOX_n == 1 ~ "1",
    EBOX_n == 2 ~ "2",
    EBOX_n >= 3 ~ ">3"
  ))

# Step 2: Define groups and EBOX categories
groups <- c("Unique", "Correlation", "All")
ebox_categories <- c("0", "1", "2", ">3")

# Step 3: Function to compute Fisher test between two groups for an EBOX category
fisher_compare_ebox <- function(group1, group2, cat) {
  data_sub <- df_stats_ebox %>%
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
ebox_results <- list()

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    for (cat in ebox_categories) {
      res <- fisher_compare_ebox(g1, g2, cat)
      ebox_results[[length(ebox_results) + 1]] <- res
    }
  }
}

# Step 5: Combine and print the result table
ebox_all_results <- do.call(rbind, ebox_results)
print(ebox_all_results)

### full knockdown statistics
library(dplyr)
library(tidyr)

# Step 1: Categorize avg_knockdown and group transcripts
df_stats_knockdown <- transcript_annotation %>%
  filter(!is.nan(avg_knockdown), avg_knockdown > 0) %>%
  mutate(
    group = case_when(
      transcript_id %in% transcript_unique$transcript_id ~ "Unique",
      transcript_id %in% transcript_correlation_all$transcript_id ~ "Correlation",
      TRUE ~ "All"
    ),
    knockdown_bin = case_when(
      avg_knockdown <= 0.25 ~ "<=0.25",
      avg_knockdown <= 0.5 ~ "0.26–0.5",
      avg_knockdown <= 0.75 ~ "0.51–0.75",
      avg_knockdown <= 1 ~ "0.76–1",
      avg_knockdown > 1 ~ ">1"
    )
  )

# Step 2: Define groups and knockdown bins
groups <- c("Unique", "Correlation", "All")
knockdown_bins <- c("<=0.25", "0.26–0.5", "0.51–0.75", "0.76–1", ">1")

# Step 3: Function for Fisher test on knockdown bins
fisher_compare_knockdown <- function(group1, group2, bin) {
  data_sub <- df_stats_knockdown %>%
    filter(group %in% c(group1, group2)) %>%
    mutate(in_bin = knockdown_bin == bin)
  
  tab <- table(data_sub$group, data_sub$in_bin)
  
  test <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  
  data.frame(
    Comparison = paste(group1, "vs", group2),
    Knockdown_Bin = bin,
    P_value = signif(test$p.value, 4),
    stringsAsFactors = FALSE
  )
}

# Step 4: Loop through all combinations
knockdown_results <- list()

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    for (bin in knockdown_bins) {
      res <- fisher_compare_knockdown(g1, g2, bin)
      knockdown_results[[length(knockdown_results) + 1]] <- res
    }
  }
}

# Step 5: Combine and print results
knockdown_all_results <- do.call(rbind, knockdown_results)
print(knockdown_all_results)

## average knockdown statistics
library(ggpubr)

# Global Kruskal-Wallis test across all groups
kruskal_test <- compare_means(avg_knockdown ~ group, data = knockdown_data, method = "kruskal.test")
print(kruskal_test)

# Pairwise Wilcoxon test with p-value adjustment
pairwise_tests <- compare_means(
  avg_knockdown ~ group,
  data = knockdown_data,
  method = "wilcox.test",
  p.adjust.method = "fdr"
)
print(pairwise_tests)

ggplot(knockdown_data, aes(x = group, y = avg_knockdown, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(method = "kruskal.test", label.y = 2.1) +  # global p
  stat_compare_means(
    comparisons = list(
      c("All", "Correlation"),
      c("All", "Unique"),
      c("Correlation", "Unique")
    ),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = c(1.9, 2.0, 2.1)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Filtered avg_knockdown per Transcript Group (0 < x ≤ 2)",
    y = "Average Knockdown"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, margin = margin(r = 15))
  )

