annotated <- readRDS("Data/hg38_annotated_transcripts_for_analysis_FINAL.rds")
final_list <- read.csv("Data/CoCor_final_list_protein_coding.csv")
tsoi_ccle_list <- read.csv("Data/CCLE_Tsoi_discordant_protein_coding.csv")

final_merged <- merge(final_list,annotated,by="transcript_id")
tsoi_ccle_merged <-   merge(tsoi_ccle_list,annotated,by="transcript_id")