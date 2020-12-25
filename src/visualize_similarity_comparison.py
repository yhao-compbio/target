# !/usr/bin/env python
# created by Yun Hao @MooreLab 2020
# This script visualizes comparison of feature similarity between target pairs with common functions and those without common functions.


# Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments
sim_file	<- 'data/target_similarity/descriptor_all_select_features_mc_0.85_function_similarity_group.tsv'
plot_file	<- 'plot/target_similarity/descriptor_all_select_features_mc_0.85_function_similarity_group_boxplot.pdf'

## 1. Obtain computed feature similarity and function groups of target pairs 
# read in data frame of model feature similarity between target pairs  
group_sim_df = pd.read_csv(sim_file, sep = '\t', header = 0)
# obtain name of annotation sources 
annot_source = group_sim_df.annotation_source.unique()
annot_name = []
for aso in annot_source:
	annot_name.append('\n'.join(aso.split(' ')))

## 2. Specify figure and font size of boxplot plot 
plt.figure(figsize = (6, 6))
plt.rc('font', size = 20)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)

## 3. Make boxplot comparing the feature similarity of target pairs with common functions, to the feature similarity of those without common functions 
ax = sns.boxplot(x = 'feature_similarity', y = 'annotation_source', hue = 'function_similarity_group', data = group_sim_df, order = annot_source, showfliers = False, notch = True, palette = 'Set3')
# show annotation sources as y ticks  
plt.yticks(range(0, len(annot_source)), annot_name)
# set range of x axis  
ax.set_xlim([-0.02, 0.4])
# add marker to show sources with significant p-values 
plt.plot([0, 0], [0, 1], marker = '*', color = 'r', linestyle = 'None', markersize = 15, label = 'w > w/o (FDR<0.05)')
# add labels and legend
ax.set(xlabel = 'Feature similarity')
ax.set(ylabel = None)
plt.legend(loc = 'lower right', frameon = False, bbox_to_anchor = (1.1, 0.98))

## 4. Save boxplot
plt.tight_layout()
sns.despine()
plt.savefig(plot_file)
plt.close()

