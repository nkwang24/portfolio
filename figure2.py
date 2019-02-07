# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:48:38 2018

@author: nkwang

Figure 2. Distribution of Q-Scores for identified CNVs
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set()

df = pd.read_csv(r'C:\Users\MVL\OneDrive\MVL\CNV Manuscript\data\qscores_data.csv')

num_plots = 3
fig, axes = plt.subplots(num_plots, 1, figsize=(5,8), constrained_layout=True)
fig.suptitle('')

for i in range(num_plots):
    p = sns.distplot(df.iloc[:,i].dropna(), kde=False, bins=7,
                     hist_kws={'range':(30,100),'alpha':0.8}, ax=axes[i])
    p.set_title(df.columns.values[i])
    p.set_xlabel('Q-Score Non-Diploid',fontsize=10)
    p.set_ylabel('Count',fontsize=10)
    p.tick_params(labelsize=10)

plt.show()

