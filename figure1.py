# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 17:01:53 2019

@author: nkwang

Figure 1. Clinical relevance and specificity of CNV pipeline
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set()

fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios':[1,1]},figsize=(11,6),constrained_layout=True)

#Pie chart 1
labels = ['Positive w/CNV\n4.3% (21)','NGS Positive\n46.6% (227)','Negative\n49.1% (239)']
values = np.array([21,227,239])

pie,texts = ax[0].pie(values, radius=0.85, labels=labels,wedgeprops=dict(width=0.85,edgecolor='w'),
               explode=[0.1,0.005,0.005],startangle=20,colors=['#0066cc','#79a6d2','#c18b8b'],
               labeldistance=1.12)

for text in texts:
    text.set_fontsize(12)

#Pie chart 2
labels = [['Q-Score>50\n47.1% (153)','Q-Score<=50\n52.9% (172)'],
          ['True Positive 10.5% (34)','False Positive 1.2% (4)','Not Tested\n35.4% (115)',
           'Not Tested 49.2% (160)','False Positive 3.1% (10)','True Positive 0.6% (2)']]
values = np.array([[34,4,115], [160,10,2]])

inner,texts = ax[1].pie(values.sum(axis=1), radius=0.7, labels=labels[0],
                      wedgeprops=dict(width=0.7, edgecolor='w'),startangle=30, explode=[0.005,0.005],
                      labeldistance=0.23, colors=['#0066cc','#b81414'])
for text in texts:
    text.set_color('white')
    text.set_fontsize(12)

outer,texts = ax[1].pie(values.flatten(), radius=0.96, labels=labels[1],
                      wedgeprops=dict(width=0.25, edgecolor='w'),startangle=30,
                      explode=[0.005]*6,labeldistance=1.07,
                      colors=['#1980e6','#80bfff','#a6bfd9','#d2acac','#f28c8c','#e61919'])

ax[0].text(0,1,'a',transform=ax[0].transAxes,fontsize=16,fontweight='bold',va='top',ha='right')
ax[1].text(0,1,'b',transform=ax[1].transAxes,fontsize=16,fontweight='bold',va='top',ha='right')

plt.show()
