# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 13:44:23 2018

@author: nkwang

Figure 4. Size and coverage depth distribution of identified
CNVs with Q-Score>50
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker
import numpy as np

# Set seaborn style
sns.set()

# Import data
df_cnv_counts = pd.read_csv(r'C:\Users\MVL\OneDrive\MVL\CNV Manuscript\data\filtered\cnv_counts_Q50.csv')
df_target_coverage = pd.read_csv(r'C:\Users\MVL\OneDrive\MVL\CNV Manuscript\data\filtered\target_coverages.csv')

# Calculate average coverage of each CNV
for i in df_cnv_counts.index:
    target_index = df_cnv_counts.at[i, 'targets'].split('..')
    cnv_size = 0
    cnv_coverage = 0
    for j in range(int(target_index[0])-1,int(target_index[1])):
        average_coverage = df_target_coverage.at[j,'average_coverage']
        target_size = df_target_coverage.at[j,'target_size']
        cnv_coverage += average_coverage * target_size
        cnv_size += target_size
    df_cnv_counts.at[i, 'average_coverage'] = cnv_coverage / cnv_size

# Plot CNVs with marginal histograms
log_bins = np.logspace(np.log(0.035),np.log(1500),50)
g = sns.JointGrid(x='CNV_size_kb', y='average_coverage', data=df_cnv_counts,
                  height=8)
g.plot_joint(sns.scatterplot, size=df_cnv_counts['Count'], hue=df_cnv_counts['Confirmation'],
             style=df_cnv_counts['Confirmation'], sizes=(70,300),
             palette={'Confirmed':'g','Unconfirmed':'b','False Positive':'orangered','False Negative':'k',' ':'k','  ':'k'},
             markers={'Confirmed':'s','Unconfirmed':'o','False Positive':'X','False Negative':'^',' ':'s','  ':'s'}, legend='full')
g.ax_marg_x.hist(df_cnv_counts['CNV_size_kb'], bins=log_bins)
g.ax_marg_y.hist(df_cnv_counts['average_coverage'], orientation='horizontal',
                 bins=np.arange(0,3000,200))

# Set title
g.fig.set_size_inches(8,6)
g.fig.suptitle('')
plt.subplots_adjust(left=0.105)
plt.subplots_adjust(bottom=0.1)

# Set axes
p = g.ax_joint
p.set_xscale('log')
p.set(xlim=(0.035,1200), xticks=[0.1,1,10,100,1000])
p.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
p.set(ylim=(0,2700))
p.set_xlabel('CNV Size in kb (log)')
p.set_ylabel('Mean Coverage Depth')
lgnd = p.legend(ncol=2, fontsize=10.1, labelspacing=0.31)
lgnd.legendHandles[5]._sizes = [0]
lgnd.legendHandles[6]._sizes = [0]
lgnd.legendHandles[8]._sizes = [40]
lgnd.legendHandles[9]._sizes = [70]
lgnd.legendHandles[10]._sizes = [100]
lgnd.legendHandles[11]._sizes = [130]
lgnd.legendHandles[12]._sizes = [165]
lgnd.legendHandles[13]._sizes = [200]

# Plot average coverage of panel
#p.axhspan(505, 739, color='k', alpha=0.2)
#p.axhspan(621, 623, color='k', alpha=0.4)

# Plot minimum criteria line
#x=np.linspace(0.035,219,1000)
#y=(-np.log(x)+5.39)*45
#plt.plot(x,y, color='r')
#plt.fill_between(x,0,y,facecolor='r',alpha=0.2)

