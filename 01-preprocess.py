import pandas as pd
import re

## read data
f="data/RNASeq Data 05-17-19.xlsx"
f="data/AllSample_Steroid_timpPoints_RPKM_Annot_RatioCalcluation.xlsx"
d1 = pd.read_excel(f, sheet_name=0, header=[0,1])
d2 = pd.read_excel(f, sheet_name=1, header=[0,1])
d1.columns = ['test_id' if col[0] == d1.columns[0][0] else re.sub(r'^(\d+)',r'\1_', "_".join(col)) for col in d1.columns]
d2.columns = ['test_id' if col[0] == d2.columns[0][0] else re.sub(r'^(\d+)',r'\1_', "_".join(col)) for col in d2.columns]
df = pd.merge(d1,d2, on='test_id')
df.columns = df.columns.str.replace('test_id', 'gene')
df.to_csv("results/rpkm_table.tsv",sep="\t")

#long_df = pd.melt(df, id_vars=['gene'], var_name='subject_info', value_name='expression')
#long_df[['subject_id', 'group', 'time']] = long_df['subject_info'].str.split('_', expand=True)
#long_df = long_df.drop(columns=['subject_info'])

