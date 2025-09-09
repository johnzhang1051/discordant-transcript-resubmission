transcript_annotation <- readRDS("hg38_annotated_transcripts_EXPANDING.rds")

transcript_unique <- read.csv("Data/Final_CCLE_TCGA_overlap.csv")
transcript_correlation_all <-read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")



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
  select(group, MITF_peak_n_cat, n, proportion)

# Reorder factor levels for plotting
mitf_distribution$MITF_peak_n_cat <- factor(mitf_distribution$MITF_peak_n_cat,
                                            levels = c("0", "1", "≥2"))


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
