# !/usr/bin/env python
# created by Yun Hao @MooreLab 2020
# This script visualizes relationship map between adverse events and predictive target features


# Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


## 0. Input arguments
data_file	= 'data/target_map/descriptor_all_all_adverse_event_select_features_data.tsv'
plot_file	= 'plot/target_map/descriptor_all_all_adverse_event_select_features_data_cluster.pdf'

## 1. Read adverse event-target dataset 
ae_tar_df = pd.read_csv(data_file, sep = '\t', header = 0)

## 2. Specify figure and font size of boxplot plot 
plt.rc('font', size = 15)
plt.rc('axes', titlesize = 15)
plt.rc('axes', labelsize = 15)
plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)
plt.rc('legend', fontsize = 15)

## 3. Make clustermap showing the adverse event ~ predictive target relationships
cg = sns.clustermap(ae_tar_df,
	method = 'average',
	metric = 'jaccard',
	figsize = (36, 12),
	cmap = 'Blues', 
	linewidths = 0.2,
	linecolor = 'Black',
	dendrogram_ratio = 0.05
	)
cg.ax_col_dendrogram.set_visible(False)
cg.cax.remove()

## 4. Save boxplot
plt.savefig(plot_file, bbox_inches = 'tight')
plt.close()

