final_protein_coding <- read.csv("NEW_discordant_protein_coding.csv")
correlatin_coding_only <- read.csv("NEW_correlation_protein_coding.csv")
transcript_to_gene <- read.csv("Data/transcript_to_gene.csv")
transcript_to_gene$transcript_id <- sub("\\..*$", "", transcript_to_gene$transcript_id)

final_protein_GENE <- merge(final_protein_coding, transcript_to_gene, by='transcript_id')
final_protein_GENE <- as.data.frame(final_protein_GENE[,-(1:3)])

correlation_protein_GENE <- merge(correlatin_coding_only, transcript_to_gene, by='transcript_id')
correlation_protein_GENE <- as.data.frame(correlation_protein_GENE[,-(1:3)])

write.csv(final_protein_GENE,"final_protein_coding_GENES")
write.csv(correlation_protein_GENE,"correlation_protein_coding_GENES")

all_protein_coding <- read.csv("Data/protein_coding_FINAL.csv")
all_protein_GENE <- merge(all_protein_coding, transcript_to_gene, by='transcript_id')
write.csv(all_protein_GENE,"all_protein_coding_GENES.csv")

write.csv(final_protein_GENE,"final_protein_coding_GENES")