import pandas as pd

## read data
f="data/RNASeq Data 05-17-19.xlsx"
f="data/AllSample_Steroid_timpPoints_RPKM_Annot_RatioCalcluation.xlsx"
import pandas as pd
d1 = pd.read_excel(f, sheet_name=0, header=[0,1])
d2 = pd.read_excel(f, sheet_name=1, header=[0,1])

d1.columns = pd.MultiIndex.from_tuples([('test_id' if col[0] == d1.columns[0][0] else col[0], col[1]) for col in d1.columns])
d2.columns = pd.MultiIndex.from_tuples([('test_id' if col[0] == d2.columns[0][0] else col[0], col[1]) for col in d2.columns])
d1.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in d1.columns.values]
d2.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in d2.columns.values]
df = pd.merge(d1,d2, on='test_id_test_id')
df.columns = df.columns.str.replace('test_id_test_id', 'gene')

## test
from statsmodels.stats.anova import AnovaRM
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols;
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def rn(df, gene):
    results = []
    gene_data = df[df['gene'] == gene].drop(columns='gene')
    # Reshape the DataFrame to long format
    long_format = gene_data.melt(var_name='sample_time', value_name='expression')
    long_format['group']=long_format["sample_time"].str.extract(r"(\D+COPD)")
    long_format['time']=long_format["sample_time"].str.extract(r"_(\d+)")
    for group in long_format['group'].unique():
        group_data = long_format[long_format['group'] == group]
        group_data_0h = group_data[group_data['time'] == "0"]
        if np.mean(group_data["expression"]) > 0.1:
            #model = ols('expression ~ C(time)', data=group_data).fit()
            #anova_table = sm.stats.anova_lm(model, typ=2)
	    tukey_results = pairwise_tukeyhsd(group_data['expression'], group_data['time'])
            i=[i for i, v in enumerate(tukey_results._results_table.data[1:]) if v[0] == "0"]
            p=np.min(tukey_results.pvalues[i])
            results.append({'gene': gene, 'group': group, 'pval': p})
    return results;

from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

gene_list = df['gene']  # Adjust the range as needed
results = []
with ThreadPoolExecutor(max_workers=10) as executor:
    results = list(tqdm(executor.map(lambda g: rn(df, g), gene_list), total=len(gene_list), desc="Processing Genes"))

df_res=pd.DataFrame([item for sublist in results for item in sublist if item])
o="out";
df_res.to_csv(f"{o}.tukey_results.tsv",sep="\t")
df.to_csv(f"{o}_table.tsv",sep="\t")


import pandas as pd
df=pd.read_csv("out_table.tsv",sep="\t")
df_sig=read_csv("out.tukey_results.tsv",sep="\t");

genes = df_sig["gene"][df_sig["pval"]<0.01].tolist()
df_filtered = df[df['gene'].isin(genes)]


df_heatmap = df_filtered.set_index('gene').T  # Transpose if necessary
df_heatmap = df_heatmap.apply(pd.to_numeric, errors='coerce')
df_heatmap.dropna(inplace=True)  # Drop rows with NaN values

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(12, 8))
sns.heatmap(df_heatmap, cmap='viridis', annot=True, fmt=".2f")
plt.title("Heatmap of Significant Genes (p < 0.01)")
plt.xlabel("Genes")
plt.ylabel("Sample Time Points")
plt.show()

results_df = pd.DataFrame([item for sublist in results for item in results])  # Flatten the list of lists
print(results_df)


results = []
with ThreadPoolExecutor(max_workers=8) as executor:  # Adjust max_workers based on your CPU
    futures = {executor.submit(rn,df,gene): gene for gene in df['gene'][1:100]}
    for future in as_completed(futures):
        gene_results = future.result()
        results.extend(gene_results)

   
import numpy as np
plt.hist(np.array(pd.DataFrame(results)["pval"].dropna()))

from scipy.stats import ttest_rel
from statsmodels.stats.multicomp import pairwise_tukeyhsdu
timepoints = ['2h', '4h', '6h']
samples = sorted({col.split('_')[0] for col in df.columns if col != 'gene'})
g=df.columns.str.extract(r"(\D+COPD)")
t=df.columns.str.extract(r"_(\d+)")

for i in ["AECOPD","StableCOPD"]:
	for j in [2,4,6]:
		
		


results = []
for gene in df['gene']:
    for sample in samples:
        # Extract 0h and other time points
        time_0h = df.loc[df['gene'] == gene, f'{sample}_0h'].values
        for tp in timepoints:
            time_tp = df.loc[df['gene'] == gene, f'{sample}_{tp}'].values
            t_stat, p_value = ttest_rel(time_0h, time_tp)
            results.append({
                'gene': gene,
                'sample': sample,
                'timepoint': tp,
                't_stat': t_stat,
                'p_value': p_value
            })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)



