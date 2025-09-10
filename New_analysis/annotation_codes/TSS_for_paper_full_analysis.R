library(tidyverse)

# Load data
final_transcript <-read.csv("Data/final_list_new.txt")
correlation_only <- read.csv("Data/correlations_list_new.txt")
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")

# Merge to get protein_coding subset
tss_data_protein_coding <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Function to process a TSS dataframe with a distance threshold
process_tss <- function(tss_df, dist_threshold) {
  tss_df %>%
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
      is_isolated = min_dist >= dist_threshold,
      group = case_when(
        transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
        transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
        TRUE ~ "all"
      )
    )
}

# Function to plot boxplot
plot_tss <- function(df, title, filename) {
  p <- ggplot(df, aes(x = group, y = min_dist, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    scale_y_continuous(limits = c(0, 2000)) +
    labs(title = title, x = "Transcript Group", y = "Min TSS Distance (bp)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  ggsave(filename, plot = p, width = 6, height = 4)
}

# Analyze all combinations
all_sets <- list(
  list(name = "All_Transcripts", data = tss_data),
  list(name = "Protein_Coding", data = tss_data_protein_coding)
)

thresholds <- c(500, 50)

for (set in all_sets) {
  for (thresh in thresholds) {
    processed <- process_tss(set$data, dist_threshold = thresh)
    plot_title <- paste(set$name, "- Threshold", thresh, "bp")
    file_name <- paste0("Plots/Boxplot_", set$name, "_", thresh, "bp.png")
    plot_tss(processed, plot_title, file_name)
  }
}


## updated
library(tidyverse)

# Load data
tss_data <- read.csv("Data/biomart_TSS_chrom_simple.csv")
transcript_type <- read.csv("Data/transcripttype.csv")
correlation_only <- read.csv("Data/correlation_only_list.txt")   # should have column transcript_id
final_transcript <- read.csv("Data/final_transcript_list.txt")   # should have column transcript_id

# Get protein coding subset
tss_data_protein_coding <- tss_data %>%
  left_join(transcript_type, by = "transcript_id") %>%
  filter(transcript_type == "protein_coding")

# Function to process TSS and calculate spacing + group
process_tss <- function(tss_df, dist_threshold) {
  tss_df %>%
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
      is_isolated = min_dist >= dist_threshold,
      group = case_when(
        transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
        transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
        TRUE ~ "all"
      )
    )
}

# Function to plot boxplot with median + sd
plot_tss <- function(df, title, filename) {
  # Summary for label placement
  stats <- df %>%
    group_by(group) %>%
    summarise(
      median_val = median(min_dist, na.rm = TRUE),
      sd_val = sd(min_dist, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(df, aes(x = group, y = min_dist, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_text(
      data = stats,
      aes(x = group, y = 480, label = paste0("Median±SD:\n", round(median_val, 1), " ± ", round(sd_val, 1))),
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_y_continuous(limits = c(0, 500)) +
    labs(title = title, x = "Transcript Group", y = "Min TSS Distance (bp)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  ggsave(filename, plot = p, width = 6, height = 4)
}

# Dataset list
all_sets <- list(
  list(name = "All_Transcripts", data = tss_data),
  list(name = "Protein_Coding", data = tss_data_protein_coding)
)

# Distance thresholds
thresholds <- c(500, 50)

# Run all combinations
for (set in all_sets) {
  for (thresh in thresholds) {
    processed <- process_tss(set$data, dist_threshold = thresh)
    plot_title <- paste(set$name, "- Threshold", thresh, "bp")
    file_name <- paste0("Plots/Boxplot_", set$name, "_", thresh, "bp.png")
    plot_tss(processed, plot_title, file_name)
  }
}


# Function to calculate isolation percentages
compute_isolation_stats <- function(df, threshold) {
  df %>%
    mutate(
      group = case_when(
        transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
        transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
        TRUE ~ "all"
      ),
      is_isolated = min_dist >= threshold
    ) %>%
    group_by(group, is_isolated) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    filter(is_isolated == TRUE) %>%
    select(group, percent) %>%
    mutate(threshold = threshold)
}

# Run this for each dataset and each threshold
# Reuse process_tss with threshold = 0 to just get min_dist
all_results <- list()

for (set in all_sets) {
  processed <- process_tss(set$data, dist_threshold = 0)  # retain all min_dist values
  
  for (thresh in c(50, 500)) {
    stats <- compute_isolation_stats(processed, threshold = thresh) %>%
      mutate(dataset = set$name)
    
    all_results[[length(all_results) + 1]] <- stats
  }
}

# Combine into one dataframe
isolation_df <- bind_rows(all_results)

# Define custom colors
custom_colors <- c(
  "all" = "#FDDDE6",
  "correlation_only" = "#F08080",
  "final_transcript" = "#C00000"
)

# Plot with custom fill
p <- ggplot(isolation_df, aes(x = group, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  facet_grid(threshold ~ dataset) +
  labs(
    title = "Percentage of Isolated Transcripts",
    y = "Percent with min_dist ≥ threshold",
    x = "Transcript Group"
  ) +
  theme_minimal()

print(p)

# --- Simplified: Statistical comparison only final_transcript vs all ---
compare_stats <- function(df) {
  for (dset in unique(df$dataset)) {
    for (thr in unique(df$threshold)) {
      subset_df <- df %>%
        filter(dataset == dset, threshold == thr, group %in% c("all", "final_transcript"))
      
      # Must have both groups to run test
      if (length(unique(subset_df$group)) == 2) {
        pct_vals <- split(subset_df$percent, subset_df$group)
        
        test <- wilcox.test(pct_vals$final_transcript, pct_vals$all)
        
        cat("\n", dset, "- Threshold", thr, "bp: final_transcript vs all\n")
        print(test)
      }
    }
  }
}

# Run test
compare_stats(isolation_df)




### TRUE/FALSE PLOTS
# Function to compute full TRUE/FALSE breakdown
compute_full_isolation_stats <- function(df, threshold, dataset_name) {
  df %>%
    mutate(
      group = case_when(
        transcript_id %in% final_transcript$transcript_id ~ "final_transcript",
        transcript_id %in% correlation_only$transcript_id ~ "correlation_only",
        TRUE ~ "all"
      ),
      is_isolated = min_dist >= threshold,
      dataset = dataset_name  # <- Add dataset name before grouping
    ) %>%
    group_by(dataset, group, is_isolated) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dataset, group) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    mutate(threshold = paste0("≥", threshold, " bp"))
}

# Run this for all datasets and thresholds
full_isolation_data <- list()

for (set in all_sets) {
  processed <- process_tss(set$data, dist_threshold = 0)
  for (thresh in c(50, 500)) {
    stats <- compute_full_isolation_stats(processed, thresh, dataset_name = set$name)
    full_isolation_data[[length(full_isolation_data) + 1]] <- stats
  }
}


# Combine all
full_isolation_df <- bind_rows(full_isolation_data)

# Plot TRUE vs FALSE
ggplot(full_isolation_df, aes(x = group, y = percent, fill = is_isolated)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(threshold ~ dataset) +
  labs(
    title = "Isolation Status Breakdown by Group",
    x = "Transcript Group",
    y = "Percent of Transcripts",
    fill = "Isolated (min_dist)"
  ) +
  theme_minimal()

