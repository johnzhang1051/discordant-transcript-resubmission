# Tsoi

* In this folder, we merge the Tsoi transcript-MITF correlations and the CCLE correlations, to ensure that the transcript expression correlations we're seeing are consistent among different datasets
* This increases confidence that our results aren't due to random occurences weird expression data

## This folder has the following steps:

1. In `CCLE_Tsoi_final_discordant.R`
  * We import data from Kallisto
  * Run pearson and spearman correlations for transcript-MITF (ENST00000394351), no gene-level correlations
  * Filter to correlations > 0.5
  * Inner join with discordant transcripts from CCLE (`/CoCor_analysis/CoCor_CCLE.R`)
2. In `MITF_correlated_final.R`
  * We merge the filtered Tsoi correlation data with the CCLE data
  * Then filter to protein-coding transcripts only
  * Inner join with discordant transcripts from CCLE (`/CoCor_analysis/CoCor_CCLE.R`)
3. Export final data as `CCLE_Tsoi_discordant_protein_coding.csv`, 


## Notes:
* I'm not sure what the results in this folder are used for
* For example, `final_MITF_correlated.csv` is not used anywhere else I can find
* I'm guessing that `CCLE_Tsoi_discordant_protein_coding.csv` is the final result, and is used in:
  * `siMITF/siMITF_no_filtering.R`, but the file is renamed as `discordant.csv`
