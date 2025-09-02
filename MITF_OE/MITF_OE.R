library(dplyr)
library(ggplot2)

GSE163646<- read.csv("Data/GSE_163646_OE_kallisto_transcript_TPM.csv")
PRJNA704810 <- read.csv("Data/PRJNA704810_OE_kallisto_transcript_TPM.csv")

merged_OE <- merge(GSE163646, PRJNA704810, by="transcript_id")
# Filter to exclude any row with expression ≤ 0.1 in any sample column
merged_OE_filtered <- merged_OE %>%
  filter(if_all(-transcript_id, ~ . > 0.1))


# Define sample groups
prjna_con <- c("SRR13782518", "SRR13782519")
prjna_oe  <- c("SRR13782520", "SRR13782521", "SRR13782522")
gse_oe    <- c("SRR13282354", "SRR13282355", "SRR13282356")
gse_con   <- c("SRR13282357", "SRR13282358", "SRR13282359")

# Step 1: Calculate group means using rowMeans
mean_expr_df <- merged_OE_filtered %>%
  mutate(
    PRJNA_CON_mean = rowMeans(select(., all_of(prjna_con)), na.rm = TRUE),
    PRJNA_OE_mean  = rowMeans(select(., all_of(prjna_oe)),  na.rm = TRUE),
    GSE_CON_mean   = rowMeans(select(., all_of(gse_con)),   na.rm = TRUE),
    GSE_OE_mean    = rowMeans(select(., all_of(gse_oe)),    na.rm = TRUE)
  ) %>%
  select(transcript_id, PRJNA_CON_mean, PRJNA_OE_mean, GSE_CON_mean, GSE_OE_mean)

# Preview
head(mean_expr_df)


# Step 2: Calculate OE/CON ratios for each dataset
ratio_df <- mean_expr_df %>%
  mutate(
    PRJNA_OE_CON_ratio = PRJNA_OE_mean / PRJNA_CON_mean,
    GSE_OE_CON_ratio   = GSE_OE_mean / GSE_CON_mean
  ) %>%
  select(transcript_id, PRJNA_OE_CON_ratio, GSE_OE_CON_ratio)

# Preview
head(ratio_df)

# Remove rows with Inf, -Inf, or NaN in either ratio column
ratio_df_clean <- ratio_df %>%
  filter(is.finite(PRJNA_OE_CON_ratio), is.finite(GSE_OE_CON_ratio))

# Preview cleaned result
head(ratio_df_clean)

# Read in the group files
protein_coding <- read.csv("Data/protein_coding_RESUBMISSION.csv")
correlated     <- read.csv("Data/correlated_RESUBMISSION.csv")
discordant     <- read.csv("Data/discordant_RESUBMISSION.csv")
library(dplyr)
correlated <- correlated %>%
  filter(!(transcript_id %in% discordant$transcript_id))

# Set mutually exclusive group assignments
ratio_df_annotated <- ratio_df_clean %>%
  mutate(group = case_when(
    transcript_id %in% correlated$transcript_id     ~ "correlated",
    transcript_id %in% discordant$transcript_id     ~ "discordant",
    transcript_id %in% protein_coding$transcript_id ~ "protein_coding",
    TRUE                                             ~ "other"
  ))


# PRJNA704610
library(ggplot2)
ggplot(ratio_df_annotated, aes(x = group, y = PRJNA_OE_CON_ratio)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "MITF_OE / CON ratio (PRJNA704610)", x = "Transcript Group", y = "Expression Ratio") +
  theme_minimal()

# GSE16346
ggplot(ratio_df_annotated, aes(x = group, y = GSE_OE_CON_ratio)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "MITF_OE / CON ratio (GSE16346)", x = "Transcript Group", y = "Expression Ratio") +
  theme_minimal()


# Step 1: Filter mean_expr_df to just discordant transcripts
discordant_means <- mean_expr_df %>%
  filter(transcript_id %in% discordant$transcript_id) %>%
  select(transcript_id, PRJNA_CON_mean, PRJNA_OE_mean, GSE_CON_mean, GSE_OE_mean)

# Step 2: View the result
head(discordant_means)

# Define transcript ID sets
discordant_ids <- discordant$transcript_id
correlated_ids <- setdiff(correlated$transcript_id, discordant_ids)
protein_ids    <- setdiff(protein_coding$transcript_id, union(discordant_ids, correlated_ids))


ratio_df_grouped <- ratio_df_clean %>%
  mutate(group = case_when(
    transcript_id %in% discordant_ids  ~ "discordant",
    transcript_id %in% correlated_ids  ~ "correlated",
    transcript_id %in% protein_ids     ~ "protein_coding",
    TRUE                               ~ "other"
  )) %>%
  filter(group %in% c("protein_coding", "correlated", "discordant"))  # Keep only the 3 groups


ratio_df_binned <- ratio_df_grouped %>%
  mutate(
    PRJNA_bin = case_when(
      PRJNA_OE_CON_ratio > 2 ~ ">2",
      PRJNA_OE_CON_ratio > 1 ~ "1–2",
      TRUE                   ~ "<=1"
    ),
    GSE_bin = case_when(
      GSE_OE_CON_ratio > 2 ~ ">2",
      GSE_OE_CON_ratio > 1 ~ "1–2",
      TRUE                 ~ "<=1"
    )
  )

library(ggplot2)

# PRJNA: count & plot
ratio_df_binned %>%
  count(group, PRJNA_bin) %>%
  ggplot(aes(x = group, y = n, fill = PRJNA_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MITF_OE / CON Ratio Bins (PRJNA704610)",
    x = "Transcript Group",
    y = "Transcript Count",
    fill = "Ratio Bin"
  ) +
  theme_minimal()

# GSE: count & plot
ratio_df_binned %>%
  count(group, GSE_bin) %>%
  ggplot(aes(x = group, y = n, fill = GSE_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MITF_OE / CON Ratio Bins (GSE16346)",
    x = "Transcript Group",
    y = "Transcript Count",
    fill = "Ratio Bin"
  ) +
  theme_minimal()

library(ggplot2)
library(dplyr)

# Count and convert to percent
prjna_percent <- ratio_df_binned %>%
  count(group, PRJNA_bin) %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100)

# Plot
ggplot(prjna_percent, aes(x = group, y = percent, fill = PRJNA_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MITF_OE / CON Ratio Bins (PRJNA704610)",
    x = "Transcript Group",
    y = "Percent of Transcripts",
    fill = "Ratio Bin"
  ) +
  theme_minimal()

# Count and convert to percent
gse_percent <- ratio_df_binned %>%
  count(group, GSE_bin) %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100)

# Plot
ggplot(gse_percent, aes(x = group, y = percent, fill = GSE_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MITF_OE / CON Ratio Bins (GSE16346)",
    x = "Transcript Group",
    y = "Percent of Transcripts",
    fill = "Ratio Bin"
  ) +
  theme_minimal()

library(ggplot2)
library(dplyr)

# Count and convert to percent
prjna_percent <- ratio_df_binned %>%
  count(group, PRJNA_bin) %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100)

# Plot
ggplot(prjna_percent, aes(x = group, y = percent, fill = PRJNA_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MITF_OE / CON Ratio Bins (PRJNA704610)",
    x = "Transcript Group",
    y = "Percent of Transcripts",
    fill = "Ratio Bin"
  ) +
  theme_minimal()



