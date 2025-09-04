# MITF_OE

* Used in Figure ...
* Takes previously identified discordant transcripts, correlated transcripts, and all protein-coding transcripts
* Checks that when MITF is overexpressed, does TPM increase?
* Using 4 different datasets:
  * GSE163646
  * PRJNA704810
* Each dataset is split into overexpressed MITF, and control (no overexpression)
* Then, we determine what the difference is with: `PRJNA_OE_mean / PRJNA_CON_mean` or `GSE_OE_mean / GSE_CON_mean`
* We then create binned boxplots
   * The bins are `>2-fold increase`, `1-2 fold increase`, or `<1 fold increase`
* Finally, for each group of transcripts (discordant, correlated, and all protein-coding), we see whether discordant transcripts increase more than the other groups
* There are 2 .R files:
   1. `MITF_OE.R` does this for transcripts where we filter to only those with >10 count in >=25% of samples
      * This logic was asked for by the paper reviewers
   2. `MITF_OE_no_filtering.R` does not apply this filtering logic


