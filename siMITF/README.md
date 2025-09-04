# siMITF

* In this folder we are performing MITF knockdown analysis
* Using 4 different datasets:
  * PRJEB30337
  * GSE163646
  * Henja
  * GSE283655
* We split each dataset into knockdown/siMITF and control
* Then, we determine what the difference is with: code such as `PRJEB30337_siMITF_mean / PRJEB30337_CON_mean`
* Finally, for each group of transcripts (discordant, correlated, and all protein-coding), we see whether discordant transcripts increase more than the other groups
* Used in figure 1
