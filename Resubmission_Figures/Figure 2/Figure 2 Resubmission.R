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
  select(group, MITF_peak_n_cat, n, proportion)

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
  select(group, EBOX_n_cat, n, proportion)
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
  select(transcript_id, group, avg_knockdown)

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

