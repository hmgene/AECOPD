# AECOPD
Bioinformatics Work - Acute Exacerbations of Chronic Obstructive Pulmonary Disease 

## Results

![Screen Shot 2025-02-18 at 5 32 35 PM](https://github.com/user-attachments/assets/7dc52ea6-d22a-4db3-86c3-154c4d0e3b40)
Fig 1. Normalized top 1k variable genes expression. 


We conducted a retest of the ANOVA analysis (p-value < 0.01), focusing on a nested model that examines the effects of time and group. The analysis revealed a significant time effect within each group, as indicated by a p-value threshold of 0.01. It is noteworthy that we applied a relaxed threshold (without multiple correction or false discovery rate adjustments) to demonstrate that the genes exhibiting significant time effects do not differentiate between the groups.

### Methodology:
- Data Acquisition: RNA-seq data was obtained from this source. (https://github.com/hmgene/AECOPD/blob/main/data/RNASeq%20Data%2005-17-19.xlsx)
- Analysis:
 - To correct for patient batch effects, gene expression at each time point was normalized to the 0h values.
 - Log2 fold changes for the comparisons 2h/0h, 4h/0h, and 6h/0h were used as input data.
 - The data were analyzed using repeated measures ANOVA (using the lme4 library in R) to assess the effects of group, time, and their interaction.
 - A multiple correction for false discovery rate (FDR < 0.05) was applied, resulting in the identification of 788 significant genes.
 - The expression patterns of these 788 genes were visualized by averaging the patient data across time points and groups.
 - The top annotation displays group patterns as lines, representing the average expression of each gene across the respective clusters.

### Conclusion:
Our analysis indicates that the 385 genes with significant time effects do not differentiate between the AECOPD and Stable COPD groups. The observed deviance does not exceed the level of patient heterogeneity.

