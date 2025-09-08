library(dplyr)
Kenny <- read.csv("Data/Kenny.csv")
Laurette <- read.csv("Data/Laurette.csv")
Webster <- read.csv("Data/Webster_promoter.csv")
Louph <- read.csv("Data/Louphrasitthiphol.csv")
# One-time fix: clean column name and remove version suffix
standardize_transcripts <- function(df) {
  df %>%
    rename(transcript_id = transcriptId) %>%
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))
}

Kenny <- standardize_transcripts(Kenny)
Laurette <- standardize_transcripts(Laurette)
Louph <- standardize_transcripts(Louph)



                    
# Define a function to filter and select columns
library(dplyr)
filter_promoters <- function(df) {
  df %>%
    filter(annotation == "Promoter (<=1kb)") %>%
    select(transcript_id, annotation)
}

# Apply to each dataset
Kenny_promoters <- filter_promoters(Kenny)
Laurette_promoters <- filter_promoters(Laurette)
Louph_promoters <- filter_promoters(Louph)
library(dplyr)
library(stringr)

# Define a function to clean up each dataset
clean_transcripts <- function(df) {
  df %>%
                            # Rename column
    mutate(transcript_id = str_replace(transcript_id, "\\.\\d+$", ""))  # Remove .1 or .23 at the end
}

# Apply to each promoter data frame
Kenny_promoters <- clean_transcripts(Kenny_promoters)
Laurette_promoters <- clean_transcripts(Laurette_promoters)
Louph_promoters <- clean_transcripts(Louph_promoters)





install.packages("VennDiagram")
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    Kenny = Kenny_promoters$transcript_id,
    Laurette = Laurette_promoters$transcript_id,
    Louph = Louph_promoters$transcript_id
  ),
  category.names = c("Kenny", "Laurette", "Louph"),
  filename = NULL,
  output = TRUE,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)


# To plot to screen
grid::grid.newpage()
grid::grid.draw(venn.plot)


protein_coding <- read.csv("Data/protein_coding.csv")
transcript_unique <- read.csv("Data/discordant.csv")
transcript_correlation_all <-read.csv("Data/correlated.csv")
transcript_correlation_all <- transcript_correlation_all %>%
  filter(!(transcript_id %in% transcript_unique$transcript_id))

library(dplyr)

# Filter out transcript_ids that are in transcript_unique
transcript_correlation_all <- transcript_correlation_all %>%
  filter(!(transcript_id %in% transcript_unique$transcript_id))



###
# Create a function to calculate proportions
calc_peak_proportion <- function(peaks_df, transcript_lists) {
  transcript_with_peaks <- unique(peaks_df$transcript_id)
  
  purrr::map_dfr(names(transcript_lists), function(group) {
    group_ids <- transcript_lists[[group]]
    n_total <- length(group_ids)
    n_overlap <- sum(group_ids %in% transcript_with_peaks)
    
    tibble(
      group = group,
      with_peak = n_overlap,
      total = n_total,
      proportion = n_overlap / n_total
    )
  })
}


###
# Create named list of transcript sets
transcript_sets <- list(
  protein_coding = protein_coding$transcript_id,
  unique = transcript_unique$transcript_id,
  correlation = transcript_correlation_all$transcript_id
)


###
# Apply to each ChIP-seq dataset
kenny_stats <- calc_peak_proportion(Kenny_promoters, transcript_sets) %>% mutate(dataset = "Kenny")
laurette_stats <- calc_peak_proportion(Laurette_promoters, transcript_sets) %>% mutate(dataset = "Laurette")
louph_stats <- calc_peak_proportion(Louph_promoters, transcript_sets) %>% mutate(dataset = "Louph")

# Combine all results
all_stats <- bind_rows(kenny_stats, laurette_stats, louph_stats)

library(ggplot2)

ggplot(all_stats, aes(x = group, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Proportion of Transcripts with ChIP-seq Peaks",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14)





### peak in any data set
# Combine transcript IDs from all peak sets
all_peaks <- unique(c(
  Kenny_promoters$transcript_id,
  Laurette_promoters$transcript_id,
  Louph_promoters$transcript_id
))

# Function for single overlap summary
compute_any_peak_proportion <- function(transcript_group, label) {
  overlap_count <- sum(transcript_group$transcript_id %in% all_peaks)
  total <- nrow(transcript_group)
  
  tibble(
    group = label,
    with_peak = overlap_count,
    total = total,
    proportion = overlap_count / total
  )
}

# Run for each group
any_peak_stats <- bind_rows(
  compute_any_peak_proportion(protein_coding, "protein_coding"),
  compute_any_peak_proportion(transcript_unique, "unique"),
  compute_any_peak_proportion(transcript_correlation_all, "correlation")
)

ggplot(any_peak_stats, aes(x = group, y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  labs(
    title = "Proportion of Transcripts with a Peak in Any Dataset",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 14)


### extend window
# Generalized filtering by promoter window
filter_promoter_window <- function(df, max_kb) {
  allowed_annotations <- switch(
    as.character(max_kb),
    "1" = "Promoter (<=1kb)",
    "2" = c("Promoter (<=1kb)", "Promoter (1-2kb)"),
    "3" = c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")
  )
  df %>% filter(annotation %in% allowed_annotations)
}

# Wrapper to calculate proportions for a given window
get_all_stats_by_window <- function(max_kb, label_suffix) {
  kenny <- filter_promoter_window(Kenny, max_kb) %>%
    select(transcript_id) %>%
    mutate(dataset = "Kenny")
  
  laurette <- filter_promoter_window(Laurette, max_kb) %>%
    select(transcript_id) %>%
    mutate(dataset = "Laurette")
  
  louph <- filter_promoter_window(Louph, max_kb) %>%
    select(transcript_id) %>%
    mutate(dataset = "Louph")
  
  # Compute stats
  transcript_sets <- list(
    protein_coding = protein_coding$transcript_id,
    unique = transcript_unique$transcript_id,
    correlation = transcript_correlation_all$transcript_id
  )
  
  bind_rows(
    calc_peak_proportion(kenny, transcript_sets) %>% mutate(dataset = "Kenny", window = paste0("<= ", max_kb, "kb")),
    calc_peak_proportion(laurette, transcript_sets) %>% mutate(dataset = "Laurette", window = paste0("<= ", max_kb, "kb")),
    calc_peak_proportion(louph, transcript_sets) %>% mutate(dataset = "Louph", window = paste0("<= ", max_kb, "kb"))
  )
}

# Combine all windows
all_window_stats <- bind_rows(
  get_all_stats_by_window(1, "1kb"),
  get_all_stats_by_window(2, "2kb"),
  get_all_stats_by_window(3, "3kb")
)

ggplot(all_window_stats, aes(x = group, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~window) +
  labs(
    title = "Proportion of Transcripts with ChIP-seq Peaks by Promoter Distance",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14)

# Function to get transcript_ids with peaks in any dataset, by window
get_any_peak_ids <- function(max_kb) {
  unique(c(
    filter_promoter_window(Kenny, max_kb)$transcript_id,
    filter_promoter_window(Laurette, max_kb)$transcript_id,
    filter_promoter_window(Louph, max_kb)$transcript_id
  ))
}

# Compute and combine for 1kb, 2kb, 3kb
any_peak_stats_all <- bind_rows(
  get_any_peak_ids(1) %>% {tibble(transcript_id = .)} %>%
    {bind_rows(
      compute_any_peak_proportion(protein_coding, "protein_coding"),
      compute_any_peak_proportion(transcript_unique, "unique"),
      compute_any_peak_proportion(transcript_correlation_all, "correlation")
    )} %>% mutate(window = "<= 1kb"),
  
  get_any_peak_ids(2) %>% {tibble(transcript_id = .)} %>%
    {bind_rows(
      compute_any_peak_proportion(protein_coding, "protein_coding"),
      compute_any_peak_proportion(transcript_unique, "unique"),
      compute_any_peak_proportion(transcript_correlation_all, "correlation")
    )} %>% mutate(window = "<= 2kb"),
  
  get_any_peak_ids(3) %>% {tibble(transcript_id = .)} %>%
    {bind_rows(
      compute_any_peak_proportion(protein_coding, "protein_coding"),
      compute_any_peak_proportion(transcript_unique, "unique"),
      compute_any_peak_proportion(transcript_correlation_all, "correlation")
    )} %>% mutate(window = "<= 3kb")
)

# Final plot: any-dataset, any-window
ggplot(any_peak_stats_all, aes(x = group, y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  facet_wrap(~window) +
  labs(
    title = "Proportion of Transcripts with a Peak in Any Dataset",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 14)



### plot peak within first 3 kb
# Filter by <=3kb window
kenny_3kb <- filter_promoter_window(Kenny, 3) %>%
  select(transcript_id)
laurette_3kb <- filter_promoter_window(Laurette, 3) %>%
  select(transcript_id)
louph_3kb <- filter_promoter_window(Louph, 3) %>%
  select(transcript_id)

# Define transcript groups again
transcript_sets <- list(
  protein_coding = protein_coding$transcript_id,
  unique = transcript_unique$transcript_id,
  correlation = transcript_correlation_all$transcript_id
)

# Compute proportions per dataset
three_kb_combined_stats <- bind_rows(
  calc_peak_proportion(kenny_3kb, transcript_sets) %>% mutate(dataset = "Kenny"),
  calc_peak_proportion(laurette_3kb, transcript_sets) %>% mutate(dataset = "Laurette"),
  calc_peak_proportion(louph_3kb, transcript_sets) %>% mutate(dataset = "Louph")
)

ggplot(three_kb_combined_stats, aes(x = group, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Proportion of Transcripts with Peak ≤3kb (By Dataset)",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 14)

# Combine all transcript IDs with peak within 3kb
any_3kb_transcripts <- unique(c(
  kenny_3kb$transcript_id,
  laurette_3kb$transcript_id,
  louph_3kb$transcript_id
))

any_3kb_stats <- tibble(
  group = c("protein_coding", "unique", "correlation"),
  total = c(
    length(protein_coding$transcript_id),
    length(transcript_unique$transcript_id),
    length(transcript_correlation_all$transcript_id)
  ),
  with_peak = c(
    sum(protein_coding$transcript_id %in% any_3kb_transcripts),
    sum(transcript_unique$transcript_id %in% any_3kb_transcripts),
    sum(transcript_correlation_all$transcript_id %in% any_3kb_transcripts)
  )
) %>%
  mutate(
    proportion = with_peak / total,
    dataset = "Any"
  )

ggplot(any_3kb_stats, aes(x = group, y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  labs(
    title = "Proportion of Transcripts with Peak ≤3kb (Any Dataset)",
    x = "Transcript Group",
    y = "Proportion with Peak"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 14)

### top hits
# Ensure all transcript_id columns are cleaned and standardized
kenny_3kb_ids <- filter_promoter_window(Kenny, 3)$transcript_id
laurette_3kb_ids <- filter_promoter_window(Laurette, 3)$transcript_id
louph_3kb_ids <- filter_promoter_window(Louph, 3)$transcript_id
shared_in_all <- Reduce(intersect, list(kenny_3kb_ids, laurette_3kb_ids, louph_3kb_ids))
unique_shared_in_all <- transcript_unique$transcript_id[
  transcript_unique$transcript_id %in% shared_in_all
]
print(unique_shared_in_all)

# Optional: save to CSV
write.csv(unique_shared_in_all, "unique_transcripts_with_peak_in_all_three.csv", row.names = FALSE)


###any two data sets
intersect_kenny_laurette <- intersect(kenny_3kb_ids, laurette_3kb_ids)
intersect_kenny_louph     <- intersect(kenny_3kb_ids, louph_3kb_ids)
intersect_laurette_louph  <- intersect(laurette_3kb_ids, louph_3kb_ids)
overlap_any_two <- unique(c(
  intersect_kenny_laurette,
  intersect_kenny_louph,
  intersect_laurette_louph
))
any_two <- transcript_unique$transcript_id[
  transcript_unique$transcript_id %in% overlap_any_two
]


### output peaks
# Ensure promoter peak sets are clean and unique
kenny_peaks <- Kenny_promoters %>% distinct(transcript_id)
laurette_peaks <- Laurette_promoters %>% distinct(transcript_id)
louph_peaks <- Louph_promoters %>% distinct(transcript_id)

# Initialize result with transcript_unique list
peak_counts_df <- transcript_unique %>%
  mutate(
    transcript_id = sub("\\.\\d+$", "", transcript_id)
  ) %>%
  distinct(transcript_id) %>%
  mutate(
    Kenny_peaks = as.integer(transcript_id %in% kenny_peaks$transcript_id),
    Laurette_peaks = as.integer(transcript_id %in% laurette_peaks$transcript_id),
    Louph_peaks = as.integer(transcript_id %in% louph_peaks$transcript_id)
  )

# Optional: View the result
head(peak_counts_df)

# Save to CSV
write.csv(peak_counts_df, "Transcript_peak_counts_per_dataset.csv", row.names = FALSE)
