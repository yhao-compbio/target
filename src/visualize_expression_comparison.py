# !/usr/bin/env python
# created by Yun Hao @MooreLab 2020
# This script visualizes comparison of differential expression values between selected target genes and all measured genes.


# Module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import seaborn as sns


## 0. Input arguments 
exp_file	= 'data/target_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_uniprot_median_tpm_adjust_all_adverse_event_select_targets_compare_de.tsv'
pv_file		= 'data/target_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_uniprot_median_tpm_adjust_all_adverse_event_select_targets_compare_pv.tsv'
plot_file	= 'plot/target_expression/descriptor_all_all_adverse_event_select_targets_de_compare.pdf'

## 1. Obtain differential expression data to be plotted
# read in data frame of differential expression values of selected target genes & all measured genes 
group_exp_df = pd.read_csv(exp_file, sep = '\t', header = 0)
# read in data frame that contains p value comparing the differential expression between the two above groups  
group_pv_df = pd.read_csv(pv_file, sep = '\t', header = 0)
all_aes = group_pv_df.adverse_event.values
all_tissues = group_pv_df.tissue.values
ae_fdr = group_pv_df.fdr.values

## 2. Set hyperparameters for plotting 
# set boxplot hyperparameters 
ae_name = []
sig_x = []
sig_y = []
col_dict = {}
dict_i = 0
col_set = plt.rcParams['axes.prop_cycle'].by_key()['color']
for i in range(0, len(all_aes)):
	# iterate by adverse event, obtain name of addverse event to be shown on x axis  
	if all_aes[i] in ['atrial_fibrillation', 'renal_injury', 'muscle_rigidity' , 'acute_hepatic_failure', 'rash_macular']:
		ae_name.append(' '.join(all_aes[i].split('_')))
	elif all_aes[i] == 'musculoskeletal_chest_pain':
		ae_name.append('musculoskeletal\nchest pain')
	else:
		aas = all_aes[i].split('_')
		aas_len = len(aas)
		if aas_len == 1:
			ae_name.append(aas[0])
		else:
			aas_c = ' '.join(aas[0:(aas_len-1)]) + '\n' + aas[aas_len-1]
			ae_name.append(aas_c)
	# set box color for the adverse event (grouped by the matched tissue)  
	if (all_tissues[i] in col_dict.keys()) == False:
		col_dict[all_tissues[i]] = col_set[dict_i]
		dict_i = dict_i + 1
	# set x and y axis of markers indicting comparison between the two groups is significant the adverse event  
	if ae_fdr[i] < 0.05:
		sig_x.append(i)
		sig_y.append(-0.1) 

# set legend hyperparameters
patchList = []
for cd_key in col_dict:
	# iterate by tissue, obtain name of tissue to be shown along with color patch
	if cd_key == 'Whole_Blood':
		key_name = 'Blood'
	elif cd_key == 'Heart_Atrial_Appendage':
		key_name = 'Atrium'
	elif cd_key == 'Heart_Left_Ventricle':
		key_name = 'Ventricle'
	else:
		key_name = cd_key.split('_')[0]
	# obtain color patch of each tissue  
	col_patch = mpatches.Patch(color = col_dict[cd_key], label = key_name)
	patchList.append(col_patch)

## 3. Specify figure and font size of boxplot 
plt.figure(figsize = (21, 9))
plt.rc('font', size = 25)
plt.rc('axes', titlesize = 30)
plt.rc('axes', labelsize = 30)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 30)
plt.rc('legend', fontsize = 25)

## 4. Make boxplot comparing the differential expression values of selected genes, to the differential expression values of all measured genes  
ax = sns.boxplot(x = 'adverse_event', y = 'differential_expression', hue = 'gene_group', data = group_exp_df, order = all_aes, showfliers = False, notch = True, palette = 'Set3')
ax.get_legend().remove()
# show adverse event names as x ticks 
plt.xticks(range(0, len(ae_name)), ae_name, rotation = 30, ha = 'right')
# set the color of boxes (grouped by matched tissues)
for i in range(0, len(all_tissues)):
	mybox = ax.artists[2*i]
	mybox.set_facecolor(col_dict[all_tissues[i]])
	mybox = ax.artists[2*i+1]
	mybox.set_facecolor(col_dict[all_tissues[i]])

# add marker to show adverse events with significant p-values 
plt.plot(sig_x, sig_y, marker = '*', color = 'r', linestyle = 'None', markersize = 15)
# add axis labels
ax.set(xlabel = None)
ax.set(ylabel = 'Expression')
# add legends showing the color patches of all tissues 
plt.legend(handles = patchList, loc = 'upper left', frameon = False, ncol = 5)
plt.plot(-0.3, 2.2, marker = '*', color = 'r', linestyle = 'None', markersize = 15)
plt.text(-0.1, 2.15, 'Selected genes (left) > All genes (right), FDR<0.05')

## 5. Save boxplot
plt.tight_layout()
sns.despine()
plt.savefig(plot_file)
plt.close()

