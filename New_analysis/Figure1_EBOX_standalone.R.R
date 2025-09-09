transcript_annotation <- readRDS("hg38_annotated_transcripts_EXPANDING.rds")

transcript_unique <- read.csv("Data/Final_CCLE_TCGA_overlap.csv")
transcript_correlation_all <-read.csv("Data/correlation_CCLE_0.5_TCGA_0.3_FINAL_NEW.csv")

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
  select(group, EBOX_n_cat, n, proportion)

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