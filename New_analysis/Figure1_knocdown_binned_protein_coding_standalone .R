
library(dplyr)
library(ggplot2)
library(tidyr)

transcript_annotation <- readRDS("hg38_annotated_transcripts_EXPANDING.rds")
protein_coding <- read.csv("Data/ALL_protein_coding_FINAL.csv")


transcript_annotation <- transcript_annotation[
  transcript_annotation$transcript_id %in% protein_coding$transcript_id,  # or just protein_coding if a vector
]


transcript_unique <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")
transcript_correlation_all <-read.csv("Data/MITF_correlated_CCLE__Tsoi.csv")



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
  select(group, knockdown_bin, n, proportion)

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
