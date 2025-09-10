library(tidyr)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(janitor)
library(regioneR)
require(ChIPpeakAnno)
require(GenomicRanges)
require(genomation)
require(TxDb.Hsapiens.UCSC.hg18.knownGene)
require(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## import primary melanoycte RNAseq transcript level counts
transcript_annotated <- read.csv(file = "biomart_TSS_chrom_simple.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)


## for each chromosome find transcripts that have duplicated TSS,then rowbind, then merge back to original list
dup1 <- transcript_annotated %>%  filter (Chromosome == "1") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup2 <- transcript_annotated %>%  filter (Chromosome == "2") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup3 <- transcript_annotated %>%  filter (Chromosome == "3") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup4 <- transcript_annotated %>%  filter (Chromosome == "4") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup5 <- transcript_annotated %>%  filter (Chromosome == "5") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup6 <- transcript_annotated %>%  filter (Chromosome == "6") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup7 <- transcript_annotated %>%  filter (Chromosome == "7") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup8 <- transcript_annotated %>%  filter (Chromosome == "8") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup9 <- transcript_annotated %>%  filter (Chromosome == "9") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup10 <-transcript_annotated %>%  filter (Chromosome == "10") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup11 <- transcript_annotated %>%  filter (Chromosome == "11") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup12 <- transcript_annotated %>%  filter (Chromosome == "12") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup13 <- transcript_annotated %>%  filter (Chromosome == "13") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup14 <- transcript_annotated %>%  filter (Chromosome == "14") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup15 <- transcript_annotated %>%  filter (Chromosome == "15") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup16 <- transcript_annotated %>%  filter (Chromosome == "16") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup17 <- transcript_annotated %>%  filter (Chromosome == "17") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup18 <- transcript_annotated %>%  filter (Chromosome == "18") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup19 <- transcript_annotated %>%  filter (Chromosome == "19") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup20 <- transcript_annotated %>%  filter (Chromosome == "20") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup21 <- transcript_annotated %>%  filter (Chromosome == "21") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup22 <- transcript_annotated %>%  filter (Chromosome == "22") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
dup23 <- transcript_annotated %>%  filter (Chromosome == "23") %>% group_by(Transcript.start..bp.) %>% filter(n() > 1)
all_dup <- rbind (dup1, dup2, dup3, dup4, dup5, dup6, dup7, dup8, dup9, dup10, dup11, dup12, dup13, dup14, dup15, dup16, dup17, dup18, dup19, dup20, dup21, dup22, dup23)
## add column "FALSE"(i.e. not a unqiue transcript), then will make all others TRUE
all_dup$unique_TSS <- 'FALSE'

## merge to transcripts annotated and then replace NA with TRUE
Transcripts_annotated_dups <- merge(x = transcript_annotated, y = all_dup, all.x = TRUE)
Transcripts_annotated_dups2 <- Transcripts_annotated_dups %>% mutate_all(~replace(., is.na(.), "TRUE"))             

##CACGTG is palindromic so present both strand
CACGTG <- read.csv(file = "Data/pwmscan_hg38_9649_36388 CACGTG.txt", row.names = NULL, head =FALSE, sep="\t", stringsAsFactors = FALSE)
## CACATG  is palindromic with CATGTG so also each site is present on both strands
CACATG <- read.csv(file = "Data/pwmscan_hg38_22687_36934 CACATG.txt", row.names = NULL, head =FALSE, sep="\t", stringsAsFactors = FALSE)
## merge CACGTG and CACATG
CACGTG_coordinates <- CACGTG[1:3]
write.table(CACGTG_coordinates, "CACGTG.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
CACGTG_GR <- readGeneric("Data/CACGTG.txt", chr = 1, start = 2, end = 3, strand = NULL,
                            meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                            remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

all_EBOXES <- rbind(CACATG,CACGTG)
all_EBOXES_cooridinates <- all_EBOXES[,1:3]
write.table(all_EBOXES_cooridinates, "EBoxes.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
Eboxes_GR <- readGeneric("Data/EBoxes.txt", chr = 1, start = 2, end = 3, strand = NULL,
                            meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                            remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MITF_peaks_GR[1:3]

## could add 2000 to the transcription start site and then merge to get an interval https://gt.rstudio.com/reference/cols_merge_range.html
## then need to match the ebox genomic coordiante to this interval
  
  
## find distance between eboxes and promoters
## https://gregglab.neuro.utah.edu/2013/01/24/retrieve-snps-in-promoters-and-flanking-sequence-for-many-genes/
## https://rpubs.com/josephuses626/findInterval
## https://cran.r-project.org/web/packages/intervals/intervals.pdf
  
## Chip_seq peak locations
MITF_peaks <- read.csv(file = "Data/MITF_CHIP_SEQ_hg38.txt", row.names = NULL, head =FALSE, sep="\t", stringsAsFactors = FALSE)

##
MITF_peaks_GR<- readGeneric("Data/MITF_CHIP_SEQ_hg38.txt", chr = 1, start = 2, end = 3, strand = NULL,
            meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
            remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MITF_peaks_GR[1:3]
## https://stackoverflow.com/questions/58117263/matching-intervals-with-values-in-another-table-in-r
## http://cbdm-01.zdv.uni-mainz.de/~jibnsale/teaching/chip-seq_exercises.html#get-the-subset-of-genes-that-have-a-peak-in-their-promoter-regions-and-write-their-gene-ids-into-an-output-file.
## able to use ebox ranges as "mini-peak" and use same tool
##https://rdrr.io/bioc/genomation/man/readGeneric.html
## https://bioconductor.github.io/BiocWorkshops/solving-common-bioinformatic-challenges-using-genomicranges.html#finding-overlaps-between-granges-objects

## https://www.biostars.org/p/268924/
## load the TXD
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## list the contents that are loaded into memory
ls('package:TxDb.Hsapiens.UCSC.hg38.knownGene')
## show the db object that is loaded by calling it's name
TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
genes <- genes(txdb)
transcripts<- transcripts(txdb)
transcripts[1:3]
##retreieving promoter regions
## example from above hits 
prom <- promoters(transcripts, upstream=2000, downstream=200)
prom[1:3]

\

## hits <- findOverlaps(query, subject, ignore.strand=TRUE)
## THIS SEEMS REALLY GOOD https://www.biostars.org/p/473162/
## https://www.biostars.org/p/46531/
## have to check of inversing the subject and query will 
## they give very different answers so will have to really figure this out
## should be able to look manually and confirm which list looks correct
## prom_reverse is what works when checking manually
MITFProm.Hits = findOverlaps(MITF_peaks_GR, prom)
MITF_prom = transcripts[subjectHits(MITFProm.Hits)]
write.table( MITF_prom, file = "MITF_prom.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

## remimport and count the number of MITF peak per transcripts
MITF_prom_2 <- read.csv(file = "Data/MITF_prom.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
MITF_transcripts_counts <- MITF_prom_2 %>% count(tx_name)
colnames(MITF_transcripts_counts)<- c("Transcript.stable.ID.version","n")

## problem is column names do not match, need to change
Transcripts_annotated_dups_MITF <- merge(x = Transcripts_annotated_dups2, y = MITF_transcripts_counts, by="Transcript.stable.ID.version", all.x = TRUE)
Transcripts_annotated_dups_MITF2 <- Transcripts_annotated_dups_MITF %>% mutate_all(~replace(., is.na(.), "0"))  

MITF_prom[1:3]
MITF_prom[1:100]

## hits <- findOverlaps(query, subject, ignore.strand=TRUE)
## if query and subject are inversed this does not work - out of bounds as query larger than subject
## seems to work but this ignores negative strand hits
EBOXProm.hits = findOverlaps(Eboxes_GR, prom)
EBOX_prom = transcripts[subjectHits(EBOXProm.hits)]
write.table(EBOX_prom, file = "EBOXPROM.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
## CACGTG only
CACGTGProm.hits = findOverlaps(CACGTG_GR, prom)
CACGTG_prom = transcripts[subjectHits(CACGTGProm.hits)]
write.table(CACGTG_prom, file = "CACGTGPROM.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
## remimport the  and count the number of CACGTG per transcript
CACGTG_prom_2 <- read.csv(file = "Data/CACGTGPROM.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
CACGTG_transcripts_counts <- CACGTG_prom_2 %>% count(tx_name)
write.table(CACGTG_transcripts_counts, file = "CACGTG_counts.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )



## remimport the  and count the number of EBOX per transcript
EBOX_prom_2 <- read.csv(file = "Data/EBOX_PROM.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
EBOX_transcripts_counts <- EBOX_prom_2 %>% count(tx_name)
colnames(EBOX_transcripts_counts)<- c("Transcript.stable.ID.version","n")
colnames(CACGTG_transcripts_counts)<- c("Transcript.stable.ID.version","n")
Transcripts_annotated_dups_MITF_EBOX <- merge(x = Transcripts_annotated_dups_MITF2, y = EBOX_transcripts_counts, by="Transcript.stable.ID.version", all.x = TRUE)
Transcripts_annotated_dups_MITF_EBOX2 <- merge(x = Transcripts_annotated_dups_MITF_EBOX, y = CACGTG_transcripts_counts, by="Transcript.stable.ID.version", all.x = TRUE)
Transcripts_annotated_dups_MITF_EBOX3 <- Transcripts_annotated_dups_MITF_EBOX2 %>% mutate_all(~replace(., is.na(.), "0"))  

colnames(Transcripts_annotated_dups_MITF_EBOX3)[11] <- "MITF_peak_n"
colnames(Transcripts_annotated_dups_MITF_EBOX3)[12] <- "EBOX_n"
colnames(Transcripts_annotated_dups_MITF_EBOX3)[13] <- "CACGTG_n"


## load siMITF data (siMITF/siCON ratio for each transcript from four published data sets)
si_MITF <- read.csv(file = "Data/siMITF_gene_transcript.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
final_annotated <- merge(x = Transcripts_annotated_dups_MITF_EBOX3, y = si_MITF, by="Transcript.stable.ID", all.x = TRUE)
saveRDS(final_annotated, file = "hg38_annotated_transcripts.rds")


  
hg38_annotated_transcripts <- readRDS("Data/hg38_annotated_transcripts.rds")

hg38_annotated_transcripts_Henja <- hg38_annotated_transcripts[,1:12]
Henja <- read.csv(file = "Data/Henja_siMITF.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
summary(Henja)
hg38_annotated_transcripts_Henja2 <- merge(x = Henja, y = hg38_annotated_transcripts_Henja, by="Transcript.stable.ID")

## below works well!!
# set up cut-off values 
breaks <- c(0,0.25,0.5,0.75,1,31021913035)
# specify interval/bin labels
tags <- c("[0-0.25)","[0.25-0.5)", "[0.5-0.75)", "[0.75-1)", "[1-31021913035)")
# bucketing values into bins
group_tags <- cut(hg38_annotated_transcripts_Henja2$Henja, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_tags)

## load Ellie lists and merge
Pearson_hits <- read.csv(file = "Data/Pearson_differential_hits.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
Pearson_annotated <- merge(x = Pearson_hits, y = final_annotated, by="Transcript.stable.ID")

Spearman_hits <- read.csv(file = "Data/Spearman_differential_hits.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
Spearman_annotated <- merge(x = Spearman_hits, y = final_annotated, by="Transcript.stable.ID")

x <- as.data.frame(table(final_annotated))


EBOX_count <- hg38_annotated_transcripts %>% count(EBOX_n)
CACGTG_count <- hg38_annotated_transcripts %>% count(CACGTG_n)
CHIP_seq_count <- hg38_annotated_transcripts %>% count(MITF_peak_n)
unique_TSS_count <- hg38_annotated_transcripts %>% count(unique_TSS)

## CACGTG count does not seem to make a big difference

EBOX_count_pearson <- Pearson_annotated %>% count(EBOX_n)
CACGTG_count_pearson <- Pearson_annotated %>% count(CACGTG_n)
CHIP_seq_count_pearson <- Pearson_annotated %>% count(MITF_peak_n)
unique_TSS_count_pearson <- Pearson_annotated %>% count(unique_TSS)
     
EBOX_count_spearman <- Spearman_annotated %>% count(EBOX_n)
CHIP_seq_count_spearman <- Spearman_annotated %>% count(MITF_peak_n)
unique_TSS_count_spearman <- Spearman_annotated %>% count(unique_TSS)  

##import CCLE pos pearson hits
CCLE_pos_pearson <- read.csv(file = "Data/CCLE_pearson_potentials_filtered.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
CCLE_pos_pearson_annotated <- merge(x = CCLE_pos_pearson, y = final_annotated, by="Transcript.stable.ID", all.x = TRUE)
EBOX_count_CCLE_pos_pearson <- CCLE_pos_pearson_annotated %>% count(EBOX_n)
CHIP_seq_count_CCLE_pos_pearson <- CCLE_pos_pearson_annotated %>% count(MITF_peak_n)
unique_TSS_count_CCLE_pos_pearson <- CCLE_pos_pearson_annotated %>% count(unique_TSS)

##import CCLE pos Spearman hits
CCLE_pos_spearman <- read.csv(file = "Data/CCLE_spearman_potentials_filtered.txt", row.names = NULL, head =TRUE, sep="\t", stringsAsFactors = FALSE)
CCLE_pos_spearman_annotated <- merge(x = CCLE_pos_spearman, y = final_annotated, by="Transcript.stable.ID", all.x = TRUE)
EBOX_count_CCLE_pos_spearman <- CCLE_pos_spearman_annotated %>% count(EBOX_n)
CHIP_seq_count_CCLE_pos_spearman <- CCLE_pos_spearman_annotated %>% count(MITF_peak_n)
unique_TSS_count_CCLE_pos_spearman <- CCLE_pos_spearman_annotated %>% count(unique_TSS)
