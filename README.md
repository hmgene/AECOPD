# AECOPD
Bioinformatics Work - Acute Exacerbations of Chronic Obstructive Pulmonary Disease 

## Results
We conducted a retest of the ANOVA analysis (p-value < 0.01), focusing on a nested model that examines the effects of time and group. The analysis revealed a significant time effect within each group, as indicated by a p-value threshold of 0.01. It is noteworthy that we applied a relaxed threshold (without multiple correction or false discovery rate adjustments) to demonstrate that the genes exhibiting significant time effects do not differentiate between the groups.

### Methodology:
- Data Acquisition: RNA-seq data was obtained from this source. (https://github.com/hmgene/AECOPD/blob/main/data/RNASeq%20Data%2005-17-19.xlsx)
- Preprocessing: To minimize batch effects related to patient variability, we normalized each time point against the baseline (0h).
- Testing: ANOVA was performed using the code available at this link.(https://github.com/hmgene/AECOPD/blob/main/02-heatmap.r)
- Summarization: Results are summarized here. (https://github.com/hmgene/AECOPD/tree/main/result)
- Visualization: Plots are available for review at this location. (https://github.com/hmgene/AECOPD/tree/main/result)

### Conclusion:
Our analysis indicates that the 385 genes with significant time effects do not differentiate between the AECOPD and Stable COPD groups. The observed deviance does not exceed the level of patient heterogeneity.

