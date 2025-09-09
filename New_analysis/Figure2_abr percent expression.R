library(readr)
CCLE_trans_expr <- read_tsv("Data/CCLE_Melanoma_transcript_expression.txt")

abr_expr <- CCLE_trans_expr[CCLE_trans_expr$gene_id == "ABR", ]
# Sum across all samples (columns 3 and onward)
abr_expr$sum_expr <- rowSums(abr_expr[, 3:ncol(abr_expr)], na.rm = TRUE)
# Total ABR expression across all transcripts
total_abr_expr <- sum(abr_expr$sum_expr)

# Add percent column
abr_expr$percent_expr <- abr_expr$sum_expr / total_abr_expr * 100
abr_expr_sorted <- abr_expr[order(-abr_expr$percent_expr), c("transcript_id", "sum_expr", "percent_expr")]
print(abr_expr_sorted)


###
EBox <- read.delim("Data/EBoxes.txt", sep = "\t", header = TRUE)
