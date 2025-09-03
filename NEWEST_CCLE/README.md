# NEWEST__CCLE

* This folder contains the R files and original datasets used to identify discordant transcripts in **CCLE data**
* Specifically, we clean data, run correlation analyses, and test significance between correlation values

## This folder has the following steps:

1. Cleans original datasets in `/Original_data_set/Clean_Newest_CCLE_Melanoma_Expression.R`
   * Saves the cleaned datasets as .csv files in `/Cleaned_data_sets`
   * Where are the original data sets from? Kallisto? Tsoi?
   * **We need URL's or sources for this data**
2. The cleaned data is copied to the folders `CoCor_Data_sets` and `Cleaned_data_sets`
3. In `/Filter_low_expressed`
   * We filter transcripts at 3 levels:
     * Filter correlated transcripts `filter_correlated.csv`
     * All protein coding transcripts `filter_protein_coding.R`
     * Discordant transcripts from CCLE `filter.R`
   * Results from these are outputted and used in eventual analyses and figures
4. Correlation analysis (pearson, spearman, CoCor, FDR) all run in `/CoCor_analysis/CoCor_CCLE.R`
   * Filter out sample `ACH-000931` because it has duplicate data
   * Calculates correlations (pearson, spearman) between transcript-MITF, gene-MITF, transcript-gene (in case dependency)
   * Uses [cocor.dep.groups.overlap](https://cran.r-project.org/web/packages/cocor/cocor.pdf) to test significance between the transcript-MITF and gene-MITF correlations
   * Created functions `run_cocor_test` and `run_cocor_test_spearman` to run on all valid transcript and gene pairs
   * Use `Benjamini & Hochberg` to control for `FDR` and returns a p-value, which we filter to below `0.05`
   * We create a final `CoCor_final_list_protein_coding.csv` filtering to:
     * Transcript-MITF correlation is >=0.5
     * Gene-MITF correlation <= 0.05
     * FDR p-value < 0.05
     * `transcript_type` = protein_coding
   * This will give us the **discordant transcripts** from the CCLE data
   * We continue to collect discordant transcripts in the `Tsoi` folder


## Notes:
* I'm assuming the other folders are just the old analysis (the 20% analysis, before we used CoCor instead)
* Where does `transcripttype.csv` come from?
* The csv's in `/Original_data_set` are all from Depmap
* `CCLE_Annotations` is an attempt to complile all data into one dataframe, believe this is for when we have are creating supplmentary figures that include all data for each transcript
