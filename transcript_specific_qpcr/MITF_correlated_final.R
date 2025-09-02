if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

library(Biostrings)

# Load sequences from all .fa files (assumes all are from same gene)
fa_files <- list.files("path_to_fasta_files", pattern = "\\.fa$", full.names = TRUE)
all_transcripts <- lapply(fa_files, readDNAStringSet)
names(all_transcripts) <- basename(fa_files)
