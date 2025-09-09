* Complete dataset for computational analysis
* All transcript features, analysis results, intermediate calculations
* For the paper, we'll include this table for just correlated and discordant transcripts

Features:
* `transcript_id`
* `model_id`
* `data_source`
* `gene_id`
* `transcript_type`
* `mitf_pearson_correlation`
* `mitf_spearman_correlation`
* `mitf_gene_pearson_correlation`
* `mitf_gene_spearman_correlation`
* `fdr`
* `counts`
* `mitf_overexpression`
* `has_chip_peak`
* `start_coordinate`
* `end_coordinate`
* `is_correlated`
* `is_discordant`
* `passes_filtering` # this is whether transcript passes `>10 count in >=25% of samples` logic
