# discordant-transcript-resubmission
Re-organized code for paper identifying MITF regulated transcript isoforms


# Analysis Flow:

1. NEWEST_CCLE:
   * Get original datasets from Depmap
   * Datasets:
      * OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected
      * OmicsExpressionTranscriptsExpectedCountProfile
      * Profile_Sample_ID
      * Melanoma_Sample_ID
      * Model - **not used anywhere?**
      * Melanoma_models - **not used anywhere?**
2. NEWEST_CCLE:
   * Clean data
   * Filter to protein coding transcripts
   * Generate correlation values (pearson + spearman) for transcripts correlated with MITF (>= 0.5)
3. Tsoi:
   * Using Tsoi datasets
   * Generate correlation values, then filter only to transcripts correlated with MITF (>= 0.5)
   * Join CCLE correlated transcripts with Tsoi correlated transcripts
4. NEWEST_CCLE:
   * Bring in the joined transcript list from before
   * Use CoCor to find only transcripts where transcript-MITF correlation showed a stastically significant stronger correlation than the gene-MITF correlation
   * This results in our list of get discordant transcripts
5. Create Supplementary Figures and Files:
   * Figure 1:
      * A subset of transcripts in melanoma have high correlation with MITF, while the parent genes of the same transcripts show lower correlation with MITF
   * Figure 2:
      * Discordantly correlated transcripts are enriched for features associated with MITF regulation and are more likely to be associated with a unique promoter
   * Supplementary Table 1:
      * Fully annotated dataset of all MITF-regulated transcript features
   * Supplementary Figure 1:
      * Outline of transcript filtering strategy
   * Supplementary Figure 2:
      * Expression of discordantly correlated genes is enriched in melanoma
   * Supplementary Figure 3:
      * Gene Ontology and pathway analysis
   * Supplementary Figure 4:
      * Discordantly correlated transcript group has lower average expression in CCLE melanoma as compared to MITF correlated transcripts
   * Supplementary File 1:
      * Complete dataset for computational analysis
      * All transcript features, analysis results, intermediate calculations
