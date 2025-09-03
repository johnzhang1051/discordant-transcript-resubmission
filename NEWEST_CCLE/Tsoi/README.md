# Tsoi

* In this folder, we merge the Tsoi transcript-MITF correlations and the CCLE correlations, to ensure that the transcript expression correlations we're seeing are consistent among different datasets
* We ran Pearson + Spearman correlations on transcript-MITF, filtered to only correlations >=0.5, then merged with the CCLE discordant transcripts
* This increases confidence that our results aren't due to random occurences weird expression data

## This folder has the following steps:

1. In `CCLE_Tsoi_final_discordant.R`
  * We import data from Kallisto
  * Run pearson and spearman correlations for transcript-MITF (ENST00000394351), no gene-level correlations
  * Filter to correlations > 0.5
  * Inner join with discordant transcripts from CCLE (`/CoCor_analysis/CoCor_CCLE.R`)
2. In `MITF_correlated_final.R`
  * We take the Tsoi transcripts with MITF-correlation >= 0.5
  * Merge it with correlations from CCLE transcripts
  * Filter to only CCLE transcripts with correlation >= 0.5
  * Filter to only protein-encoding transcripts
  * Export as `final_MITF_correlated.csv`
3. Export final data as `CCLE_Tsoi_discordant_protein_coding.csv`, 


## Notes:
* `CCLE_Tsoi_discordant_protein_coding.csv` is the final result, and is used in:
  * `siMITF/siMITF_no_filtering.R`, but the file is renamed as `discordant.csv`
