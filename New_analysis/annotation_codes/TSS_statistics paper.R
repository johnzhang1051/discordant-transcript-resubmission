library(tidyverse)  # loads ggplot2, dplyr, tidyr, etc.
library(ggplot2)    # if you only need ggplot2

# Load data
final_transcript <- read.csv("Data/final_transcript_list.txt")
correlation_only <- read.csv("Data/correlation_only_list.txt")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")

# Filter to protein_coding
tss_data_protein_coding <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Process with 500bp threshold and assign groups
processed <- tss_data_protein_coding %>%
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  ) %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),
    next_dist = abs(lead(TSS) - TSS),
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_isolated = min_dist >= 500,
    group = case_when(
      transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
      transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
      TRUE ~ "all"
    )
  ) %>%
  filter(!is.na(is_isolated))  # remove NA values

# Count TRUE/FALSE per group
count_data <- processed %>%
  filter(group %in% c("all", "correlation_only", "final_transcript")) %>%
  group_by(group, is_isolated) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

# Fisher's tests
table_all <- table(processed %>% filter(group %in% c("all", "final_transcript")) %>% select(group, is_isolated))
table_corr <- table(processed %>% filter(group %in% c("correlation_only", "final_transcript")) %>% select(group, is_isolated))

p_label_all <- ifelse((p_all <- fisher.test(table_all)$p.value) < 0.001, "p < 0.001", paste0("p = ", signif(p_all, 3)))
p_label_corr <- ifelse((p_corr <- fisher.test(table_corr)$p.value) < 0.001, "p < 0.001", paste0("p = ", signif(p_corr, 3)))

# Define color shades for each group + status
fill_colors <- c(
  "all_TRUE" = "#FDDDE6",
  "all_FALSE" = "#FAEBEF",
  "correlation_only_TRUE" = "#F08080",
  "correlation_only_FALSE" = "#F6B8B8",
  "final_transcript_TRUE" = "#C00000",
  "final_transcript_FALSE" = "#E06666"
)

# Add fill_key for mapping color
count_data <- count_data %>%
  mutate(fill_key = paste(group, ifelse(is_isolated, "TRUE", "FALSE"), sep = "_"))

# Compute annotation heights
max_y <- max(count_data$percent)
p_line_1 <- max_y + 0.1 * max_y
p_text_1 <- p_line_1 + 3
p_line_2 <- p_line_1 + 7
p_text_2 <- p_line_2 + 3

# Final plot
ggplot(count_data, aes(x = group, y = percent, fill = fill_key)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), aes(group = is_isolated), width = 0.6) +
  geom_text(
    aes(label = paste0(round(percent, 1), "%"), group = is_isolated),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_manual(values = fill_colors) +
  # Bracket: final_transcript vs all
  annotate("segment", x = 1 - 0.2, xend = 3 - 0.2, y = p_line_1, yend = p_line_1, size = 0.4) +
  annotate("text", x = 2 - 0.2, y = p_text_1, label = paste("final vs all:", p_label_all), size = 4, fontface = "italic") +
  # Bracket: final_transcript vs correlation_only
  annotate("segment", x = 2 + 0.2, xend = 3 + 0.2, y = p_line_2, yend = p_line_2, size = 0.4) +
  annotate("text", x = 2.5 + 0.2, y = p_text_2, label = paste("final vs corr:", p_label_corr), size = 4, fontface = "italic") +
  labs(
    title = "Percent Isolated (≥500bp) — Protein Coding",
    x = "Transcript Group",
    y = "Percent Isolated"
  ) +
  coord_cartesian(ylim = c(0, p_text_2 + 10)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")




### update with ggpatern
library(tidyverse)
library(ggpattern)

# Load data
final_transcript <- read.csv("Data/final_transcript_list.txt")
correlation_only <- read.csv("Data/correlation_only_list.txt")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")

# Filter to protein_coding
tss_data_protein_coding <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Process with 500bp threshold and assign groups
processed <- tss_data_protein_coding %>%
  mutate(
    TSS = as.integer(TSS),
    Chromosome = as.factor(Chromosome)
  ) %>%
  arrange(Chromosome, TSS) %>%
  group_by(Chromosome) %>%
  mutate(
    prev_dist = abs(TSS - lag(TSS)),
    next_dist = abs(lead(TSS) - TSS),
    min_dist = pmin(prev_dist, next_dist, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_isolated = min_dist >= 500,
    group = case_when(
      transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
      transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
      TRUE ~ "all"
    )
  ) %>%
  filter(!is.na(is_isolated))  # remove NA values

# Count TRUE/FALSE per group
count_data <- processed %>%
  filter(group %in% c("all", "correlation_only", "final_transcript")) %>%
  group_by(group, is_isolated) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

# Fisher's tests
table_all <- table(processed %>% filter(group %in% c("all", "final_transcript")) %>% select(group, is_isolated))
table_corr <- table(processed %>% filter(group %in% c("correlation_only", "final_transcript")) %>% select(group, is_isolated))

p_label_all <- ifelse((p_all <- fisher.test(table_all)$p.value) < 0.001, "p < 0.001", paste0("p = ", signif(p_all, 3)))
p_label_corr <- ifelse((p_corr <- fisher.test(table_corr)$p.value) < 0.001, "p < 0.001", paste0("p = ", signif(p_corr, 3)))

# Fill colors per group
group_colors <- c(
  "all" = "#FDDDE6",
  "correlation_only" = "#F08080",
  "final_transcript" = "#C00000"
)

# Add pattern column: TRUE = solid, FALSE = stripe
count_data <- count_data %>%
  mutate(
    pattern = ifelse(is_isolated, "none", "stripe"),
    group = factor(group, levels = c("all", "correlation_only", "final_transcript"))
  )

# Max y for bracket placement
max_y <- max(count_data$percent)
p_line_1 <- max_y + 0.1 * max_y
p_text_1 <- p_line_1 + 3
p_line_2 <- p_line_1 + 7
p_text_2 <- p_line_2 + 3

# Plot with patterns!
ggplot(count_data, aes(x = group, y = percent,
                       fill = group, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.6,
    color = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.04,
    pattern_key_scale_factor = 0.6
  ) +
  geom_text(
    aes(label = paste0(round(percent, 1), "%"), group = is_isolated),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  # Brackets
  annotate("segment", x = 1 - 0.2, xend = 3 - 0.2, y = p_line_1, yend = p_line_1, size = 0.4) +
  annotate("text", x = 2 - 0.2, y = p_text_1, label = paste("final vs all:", p_label_all), size = 4, fontface = "italic") +
  annotate("segment", x = 2 + 0.2, xend = 3 + 0.2, y = p_line_2, yend = p_line_2, size = 0.4) +
  annotate("text", x = 2.5 + 0.2, y = p_text_2, label = paste("final vs corr:", p_label_corr), size = 4, fontface = "italic") +
  labs(
    title = "Percent Isolated (≥500bp) — Protein Coding",
    x = "Transcript Group",
    y = "Percent Isolated"
  ) +
  coord_cartesian(ylim = c(0, p_text_2 + 10)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

