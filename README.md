# discordant-transcript-resubmission
Re-organized code for paper identifying MITF regulated transcript isoforms

# Data Sources:
Depmap:
   * OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv
   * OmicsExpressionTranscriptsExpectedCountProfile.csv
   * Profile_Sample_ID.csv
   * Melanoma_Sample_ID.csv

Tsoi:
   * gencode.v44.annotation.gtf

# Analysis Flow/Steps:

1. NEWEST_CCLE:
   * Download original datasets from Depmap and Tsoi
2. NEWEST_CCLE:
   * Transcript and gene expression data is cleaned in `NEWEST_CCLE/CoCor_Data_sets/Clean_Newest_CCLE_Melanoma_Expression.R`
* Discordant and correlated transcripts are idenitified in `NEWEST_CCLE/CoCor_CCLE.R` - where we calculate correlation values (pearson + spearman) between transcripts and `ENST00000394351`/MITF 
      * FDR values are also calculated using CoCor
      * We create a final `CoCor_final_list_protein_coding.csv` which consists of transcripts that have:
         * Transcript-MITF correlation (pearson or spearman) >=0.5
         * Gene-MITF correlation (pearson or spearman) <= 0.05
         * FDR p-value < 0.05
         * `transcript_type` = protein_coding
      * This will give us the **discordant transcripts** from the CCLE data
      * We continue identifying discordant transcripts using `Tsoi` data
3. Tsoi:
   * In `NEWEST_CCLE/CCLE_Tsoi_final_discordant.R`
   * Tsoi datasets are imported and organized
   * Tsoi correlations are calculated (pearson + spearman)
   * Transcripts are filtered to only those correlated with MITF (>= 0.5 for either pearson or spearman)
   * We find the overlap between the Tsoi correlated with the CCLE discordant transcripts with Tsoi correlated transcripts to make final list
   * This gives us our pre-filtered list of discordant transcripts
4. NEWEST_CCLE:
   * We then identify **"MITF-Correlated"** transcripts in `NEWEST_CCLE/MITF_correlated_final.R`
      * Here, we take transcripts from CCLE and Tsoi where Pearson or Spearman >= 0.5 (these correlations were calculated in the previous steps)
      * Then combine CCLE + Tsoi lists into one, and export as a pre-filtered list of correlated transcripts called `final_MITF_correlated.csv`
   * In `NEWEST_CCLE/Filter_low_expressed`, we run the filtering logic `>10 count in >=25% of samples` on 3 lists:
      * `NEWEST_CCLE/Filter_low_expressed/filter_discordant.R`
      * `NEWEST_CCLE/Filter_low_expressed/filter_correlated.R`
      * `NEWEST_CCLE/Filter_low_expressed/filter_protein_coding.R`
   * That gives us the final 3 lists of transcripts we'll use to create Figures
5. Resubmission_Figures:
   * Create Supplementary Figures and Files:
   * Figure 1:
      * A subset of transcripts in melanoma have high correlation with MITF, while the parent genes of the same transcripts show lower correlation with MITF
      * This was already re-done by Steve
   * Figure 2:
      * Discordantly correlated transcripts are enriched for features associated with MITF regulation and are more likely to be associated with a unique promoter
   * Supplementary Table 1:
      * Fully annotated dataset of all MITF-regulated transcript features
   * Supplementary Figure 1:
      * Outline of transcript filtering strategy
   * Supplementary Figure 2:
      * Expression of discordantly correlated genes is enriched in melanoma
      * This was already re-done by Steve
   * Supplementary Figure 3:
      * Gene Ontology and pathway analysis
   * Supplementary Figure 4:
      * Discordantly correlated transcript group has lower average expression in CCLE melanoma as compared to MITF correlated transcripts
   * Supplementary File 1:
      * Complete dataset for computational analysis
      * All transcript features, analysis results, intermediate calculations
