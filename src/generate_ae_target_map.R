# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script generates map between adverse event and predictive targets.


## functions
source("src/functions.R");


## 0. Input arguments 
ae_map_file	<- "data/target_expression/adverse_event_tissue_map.tsv";
perf_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_all_adverse_event_feature_select_implementation/descriptor_all_all_adverse_event_1_testing_performance_summary_auc_ci.tsv";

select_file	<- "data/target_disease/descriptor_all_all_adverse_event_select_features_symbol_all.tsv";
auc_cut		<- 0.65;
out_file	<- "data/target_map/descriptor_all_all_adverse_event_select_features";

## 1. Read in selected targets of adverse events 
# read in model performance of adverse events 
ae_perf <- read.delim(file = perf_file, sep = "\t", header = T);
# identify adverse events with performance above threshold 
perf_id <- which(ae_perf$target_select_model_auc > auc_cut);
perf_aes <- ae_perf$X[perf_id];
# read in data of selected features
ae_select_features <- read.delim(file = select_file, sep = "\t", header = T); 
# obtain the selected features of each adverse event of interest  
pa_id <- which(ae_select_features$adverse_event %in% perf_aes);
ae_select_features <- ae_select_features[pa_id, ];
ae_features <- group.vector.by.categories(ae_select_features$adverse_event, ae_select_features$select_target_symbol);

## 2. Generate adverse event-predictive target dataset  
# gather selected targets of all adverse events   
all_targets <- unique(unlist(ae_features));
# generate adverse event-predictive target data 
ae_target_mat <- mapply(function(af){
	tar_vec <- rep(0, length(all_targets));
	names(tar_vec) <- all_targets;	
	tar_vec[af] <- 1;
	return(tar_vec);
}, ae_features);
ae_target_mat <- t(ae_target_mat);
# remove targets selected for fewer than 5 adverse events  
col_id <- which(colSums(ae_target_mat) >= 5);
ae_target_mat <- ae_target_mat[, col_id];
# output data frame of adverse event-predictive target 
rownames(ae_target_mat) <- sapply(rownames(ae_target_mat), function(ratm) paste(strsplit(ratm, "_")[[1]], collapse = " "))
write.table(ae_target_mat, file = paste(out_file, "_data.tsv", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F);

