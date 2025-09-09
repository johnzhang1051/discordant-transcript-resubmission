### Figure 1A Pearsons


library(ggplot2)
library(dplyr)
Gene_transcript_pearson <- read.csv("Data/merged_pearsons_annotated.csv")
Gene_transcript_pearson$transcript_id <- sub("\\..*$", "", Gene_transcript_pearson$transcript_id)
Gene_transcript_pearson$gene_id <- sub("\\..*$", "", Gene_transcript_pearson$gene_id)



# Define the transcript IDs
blue_transcript <- "ENST00000544583"

red_transcripts <- c(
  "ENST00000291107", "ENST00000302538", "ENST00000536794", "ENST00000543210",
  "ENST00000570441", "ENST00000570525", "ENST00000570688", "ENST00000571022",
  "ENST00000571120", "ENST00000571306", "ENST00000571383", "ENST00000571543",
  "ENST00000571797", "ENST00000571945", "ENST00000572152", "ENST00000572441",
  "ENST00000572585", "ENST00000572650", "ENST00000573325", "ENST00000573559",
  "ENST00000573667", "ENST00000573895", "ENST00000574048", "ENST00000574139",
  "ENST00000574257", "ENST00000574266", "ENST00000574437", "ENST00000574544",
  "ENST00000574632", "ENST00000574875", "ENST00000575770", "ENST00000575934",
  "ENST00000576668", "ENST00000576964", "ENST00000577052"
)

# Add highlight column
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = case_when(
    transcript_id == blue_transcript ~ "Blue",
    transcript_id %in% red_transcripts ~ "Red",
    TRUE ~ "Other"
  ))

# Plot
ggplot() +
  # Background points (small and transparent)
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Other"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "gray70", size = 0.1, alpha = 0.3
  ) +
  # Highlighted points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight != "Other"),
    aes(x = pearson_gene, y = pearson_transcript, color = highlight),
    size = 2.0
  ) +
  scale_color_manual(values = c("Blue" = "red", "Red" = "blue")) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()


### highlight all transcripts
library(ggplot2)
library(dplyr)

# Load and clean IDs
Gene_transcript_pearson <- read.csv("Data/merged_pearsons_annotated.csv")
Gene_transcript_pearson$transcript_id <- sub("\\..*$", "", Gene_transcript_pearson$transcript_id)
Gene_transcript_pearson$gene_id <- sub("\\..*$", "", Gene_transcript_pearson$gene_id)

# Your list of transcripts to highlight in red
highlighted_transcripts <- c(
  "ENST00000254630", "ENST00000292114", "ENST00000296126", "ENST00000334805",
  "ENST00000337109", "ENST00000357311", "ENST00000360019", "ENST00000369102",
  "ENST00000395694", "ENST00000400288", "ENST00000410045", "ENST00000412952",
  "ENST00000413554", "ENST00000417581", "ENST00000443149", "ENST00000445189",
  "ENST00000455666", "ENST00000466565", "ENST00000469470", "ENST00000472419",
  "ENST00000473116", "ENST00000473200", "ENST00000474473", "ENST00000476195",
  "ENST00000477199", "ENST00000477335", "ENST00000481227", "ENST00000488772",
  "ENST00000492603", "ENST00000492627", "ENST00000495510", "ENST00000498454",
  "ENST00000508384", "ENST00000518371", "ENST00000519578", "ENST00000524817",
  "ENST00000534662", "ENST00000535603", "ENST00000536974", "ENST00000537643",
  "ENST00000537952", "ENST00000544583", "ENST00000545993", "ENST00000546655",
  "ENST00000550042", "ENST00000555814", "ENST00000562557", "ENST00000562961",
  "ENST00000563012", "ENST00000568432", "ENST00000570101", "ENST00000572073",
  "ENST00000582769", "ENST00000586057", "ENST00000591285", "ENST00000592761",
  "ENST00000605917", "ENST00000617539", "ENST00000620867"
)

# Annotate highlight column
Gene_transcript_pearson <- Gene_transcript_pearson %>%
  mutate(highlight = ifelse(transcript_id %in% highlighted_transcripts, "Highlighted", "Other"))

# Plot
ggplot() +
  # Background points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Other"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "gray70", size = 0.1, alpha = 0.3
  ) +
  # Highlighted points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Highlighted"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "red", size = 1.0
  ) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()

## highlight all 0.5 pearson as well
library(dplyr)
library(ggplot2)

ggplot(Gene_transcript_pearson, aes(x = pearson_gene, y = pearson_transcript)) +
  geom_point(
    color = "gray60", size = 0.1, alpha = 0.3
  ) +
  # Previously highlighted points
  geom_point(
    data = Gene_transcript_pearson %>% filter(highlight == "Highlighted"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "red", size = 1.0
  ) +
  # New medium pink highlights for pearson_transcript ≥ 0.5 but not already highlighted
  geom_point(
    data = Gene_transcript_pearson %>%
      filter(pearson_transcript >= 0.5 & highlight != "Highlighted"),
    aes(x = pearson_gene, y = pearson_transcript),
    color = "salmon", size = 0.1
  ) +
  labs(
    title = "Gene vs Transcript Pearson Correlation",
    x = "Gene Pearson Correlation",
    y = "Transcript Pearson Correlation"
  ) +
  theme_minimal()




### ABR Gene expression CCLE

library(ggplot2)
library(dplyr)
library(forcats)  # for fct_reorder
## 
CCLE_gene_for_figure <- readRDS("CCLE_gene_disease_Figure1B.rds")
# Prepare the data
abr_df <- CCLE_gene_for_figure %>%
  select(DepmapModelType = 3, ABR) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ABR, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ABR, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "blue", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )






### ABR TRANSCRIPT EXPRESSION CCLE
### import df and plot ABR
filtered_df_for_ABR_CCLE_transcript <- readRDS("Figure_1C_CCLE_transcript.rds")


# Prepare the data
abr_df <- filtered_df_for_ABR_CCLE_transcript %>%
  select(DepmapModelType = 3, ENST00000544583) %>%
  mutate(color_group = ifelse(DepmapModelType == "SKCM", "Melanoma", "Other"))

# Reorder DepmapModelType factor by median ABR (highest on right)
abr_df <- abr_df %>%
  mutate(DepmapModelType = fct_reorder(DepmapModelType, ENST00000544583, .fun = median))

# Plot
ggplot(abr_df, aes(x = DepmapModelType, y = ENST00000544583, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "DepMap Model Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )





#### ABR_gene_expression_TCGA
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
### load figure specific expression object
TCGA_gene_expr_annotated_with_disease <- readRDS("TCGA_gene_disease_Figure1B.rds")

# Vector of primary diseases to keep
keep_diseases <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA",
                   "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD",
                   "UCS", "UCEC", "HNSC", "SKCM")

# Extract disease labels from row 1
sample_disease <- as.character(TCGA_gene_expr_annotated_with_disease[1, 3:ncol(TCGA_gene_expr_annotated_with_disease)])
names(sample_disease) <- colnames(TCGA_gene_expr_annotated_with_disease)[3:ncol(TCGA_gene_expr_annotated_with_disease)]

# Subset to samples with desired primary diseases
samples_to_keep <- names(sample_disease)[sample_disease %in% keep_diseases]

# Clean expression dataframe
expr_clean <- TCGA_gene_expr_annotated_with_disease[1:54, ]
saveRDS(expr_clean, "TCGA_gene_expr_clean.rds")
# Filter for ABR
abr_expr <- expr_clean %>%
  filter(Gene == "ABR") %>%
  select(gene_id, Gene, all_of(samples_to_keep))

# Pivot to long format
abr_long <- abr_expr %>%
  pivot_longer(
    cols = -c(gene_id, Gene),
    names_to = "Sample",
    values_to = "Expression"
  )

# Add disease label and color group
abr_long <- abr_long %>%
  mutate(Disease = sample_disease[Sample],
         color_group = ifelse(Disease == "SKCM", "Melanoma", "Other"),
         Expression = as.numeric(Expression))

# Reorder diseases by median ABR expression (highest on right)
abr_long <- abr_long %>%
  mutate(Disease = fct_reorder(Disease, Expression, .fun = median))

# Plot
ggplot(abr_long, aes(x = Disease, y = Expression, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "blue", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "Tumor Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )









#### ABR_transcript_expression_TCGA
### ABR expression
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
TCGA_transcript_expr_annotated_with_disease <- readRDS("TCGA_transcript_disease_Figure1B.rds")
# Vector of primary diseases to keep
keep_diseases <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA",
                   "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD",
                   "UCS", "UCEC", "HNSC", "SKCM")

# Extract disease labels from row 1
sample_disease <- as.character(TCGA_transcript_expr_annotated_with_disease[1, 3:ncol(TCGA_transcript_expr_annotated_with_disease)])
names(sample_disease) <- colnames(TCGA_transcript_expr_annotated_with_disease)[3:ncol(TCGA_transcript_expr_annotated_with_disease)]

# Subset to samples with desired primary diseases
samples_to_keep <- names(sample_disease)[sample_disease %in% keep_diseases]

# Clean expression dataframe
expr_clean <- TCGA_transcript_expr_annotated_with_disease[1:54, ]
saveRDS(expr_clean,"TCGA_expr_clean.rds")

# Filter for ABR
abr_expr <- expr_clean %>%
  filter(Gene == "ABR") %>%
  select(transcript_id, Gene, all_of(samples_to_keep))

# Pivot to long format
abr_long <- abr_expr %>%
  pivot_longer(
    cols = -c(transcript_id, Gene),
    names_to = "Sample",
    values_to = "Expression"
  )

# Add disease label and color group
abr_long <- abr_long %>%
  mutate(Disease = sample_disease[Sample],
         color_group = ifelse(Disease == "SKCM", "Melanoma", "Other"),
         Expression = as.numeric(Expression))

# Reorder diseases by median ABR expression (highest on right)
abr_long <- abr_long %>%
  mutate(Disease = fct_reorder(Disease, Expression, .fun = median))

# Plot
ggplot(abr_long, aes(x = Disease, y = Expression, fill = color_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression Across Tumor Types",
    x = "Tumor Type (ordered by median ABR expression)",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

### update to plot all cancer types vs SKCM
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Keep same initial steps
TCGA_transcript_expr_annotated_with_disease <- readRDS("TCGA_transcript_disease_Figure1B.rds")

keep_diseases <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA",
                   "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD",
                   "UCS", "UCEC", "HNSC", "SKCM")

sample_disease <- as.character(TCGA_transcript_expr_annotated_with_disease[1, 3:ncol(TCGA_transcript_expr_annotated_with_disease)])
names(sample_disease) <- colnames(TCGA_transcript_expr_annotated_with_disease)[3:ncol(TCGA_transcript_expr_annotated_with_disease)]

samples_to_keep <- names(sample_disease)[sample_disease %in% keep_diseases]

expr_clean <- TCGA_transcript_expr_annotated_with_disease[1:54, ]

abr_expr <- expr_clean %>%
  filter(Gene == "ABR") %>%
  select(transcript_id, Gene, all_of(samples_to_keep))

abr_long <- abr_expr %>%
  pivot_longer(
    cols = -c(transcript_id, Gene),
    names_to = "Sample",
    values_to = "Expression"
  )

abr_long <- abr_long %>%
  mutate(
    Disease = sample_disease[Sample],
    Group = ifelse(Disease == "SKCM", "SKCM", "Other"),
    Expression = as.numeric(Expression)
  )

# Plot just SKCM vs Other
ggplot(abr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("SKCM" = "red", "Other" = "gray70")) +
  labs(
    title = "ABR Expression in SKCM vs. Other Cancers",
    x = "",
    y = "ABR Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

