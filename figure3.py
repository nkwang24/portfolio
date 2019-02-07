# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:48:38 2018

@author: nkwang

Figure 3. Detection of deletion controls identified by gel analysis
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set()

df = pd.read_csv(r'C:\Users\MVL\OneDrive\MVL\CNV Manuscript\data\unfiltered\DATA.same_filtered.txt', sep='\t')
df_zscore = pd.read_csv(r'C:\Users\MVL\OneDrive\MVL\CNV Manuscript\data\unfiltered\DATA.PCA_normalized.filtered.sample_zscores.txt', sep='\t')

titles = ['CLN3 (0/6 detected)','TRPM1 (3/3 detected)','CACNA2D4 (3/3 detected)',
          'OCA2 (0/5 detected)','PRPF31','EYS']
region = [(6226,6240),(5808,5835),(5074,5111),(5785,5807),(7001,7013),(3111,3152)]
sample_list = [['18-2264','17-1117','17-2519','18-0337','18-0937','18-2388'],
               ['17-1519','18-1466','18-2433'],
               ['17-1807','18-1403','18-2433'],
               ['17-1278','17-2656','17-3137','18-1203','18-2990'],
               ['18-2275','18-2785','18-1389'],
               ['18-2329','18-2711']]
CNV_region = [(6.5,8.5),(1.5,8.5),(16.5,26.5),(6.5,7.5),(0,0),(0,0)]

df.columns.values[0]='Sample'
#df.columns.values[1:]=list(range(1,len(df.columns)))
df['Sample'], df['Run'] = df['Sample'].str.split('.').str
df_zscore.columns.values[0]='Sample'
#df.columns.values[1:]=list(range(1,len(df.columns)))
df_zscore['Sample'], df_zscore['Run'] = df_zscore['Sample'].str.split('.').str

num_plots=4

fig, axes = plt.subplots(2, num_plots, figsize=(18,8.5), constrained_layout=True, sharey='row', sharex='col')
fig.suptitle('')
#plt.xticks(rotation=-15)

for i in range(num_plots):
    df['Sample CNV'] = 'Negative'
    for j in range(len(df)):
        if df.iloc[j,0] in sample_list[i]:
            df.iloc[j,-1] = 'Positive'
    df_region = df.iloc[:,[0, -2, -1]+list(range(region[i][0],region[i][1]+1))]
    if i == 3:
        df_region.columns.values[3:]=list(range(2,len(df_region.columns[3:])+2))[::-1]
    else:
        df_region.columns.values[3:]=list(range(1,len(df_region.columns[3:])+1))[::-1]
    tidy_df = pd.melt(df_region, id_vars=['Sample', 'Run', 'Sample CNV'],
                      var_name='ROI', value_name='Coverage Depth')
    p = sns.lineplot(x='ROI', y='Coverage Depth', data=tidy_df,
                        units='Sample', estimator=None,
                        hue='Sample CNV', size='Sample CNV',
                        sizes={'Negative':0.2,'Positive':1},
                        ax=axes[0,i],legend=False)
    p.set_title(titles[i], fontsize=13)
    p.axvspan(CNV_region[i][0], CNV_region[i][1], color='red', alpha=0.25)
    
#axes[0,0].set_yscale('log')
#axes[0,0].set_yticks([5,20,50,125,300,750,1500,3000])
#axes[0,0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#axes[0,0].set_ylim(bottom=10)
axes[0,0].set_ylim(top=2500)

for i in range(num_plots):
    df_zscore['Sample CNV'] = 'Negative'
    for j in range(len(df_zscore)):
        if df_zscore.iloc[j,0] in sample_list[i]:
            df_zscore.iloc[j,-1] = 'Positive'
    df_region = df_zscore.iloc[:,[0, -2, -1]+list(range(region[i][0],region[i][1]+1))]
    if i == 3:
        exons = list(range(2,len(df_region.columns[3:])+2))
    else:
        exons = list(range(1,len(df_region.columns[3:])+1))
    df_region.columns.values[3:]=exons[::-1]
    if i == 1 or i == 2:
        exon_labels = exons[0::3]
    else:
        exon_labels = exons[0::2]
    tidy_df = pd.melt(df_region, id_vars=['Sample', 'Run', 'Sample CNV'],
                      var_name='ROI', value_name='Normalized Z-Score')
    p = sns.lineplot(x='ROI', y='Normalized Z-Score', data=tidy_df,
                    units='Sample', estimator=None,
                    hue='Sample CNV', size='Sample CNV',
                    sizes={'Negative':0.2,'Positive':1},
                    ax=axes[1,i],legend=False)
    p.set_xticks(exon_labels)
    p.set_xlabel('Exons')
    p.axvspan(CNV_region[i][0], CNV_region[i][1], color='red', alpha=0.25)

axes[1,0].set_ylim(-20,15)

print("Done")