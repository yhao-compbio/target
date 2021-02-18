# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script implements pvclust package to identify significant clusters of adverse events from adverse event-predictive targets map.


## function  
library(pvclust);


## 0. Input arguments
map_file	<- "data/target_map/descriptor_all_all_adverse_event_select_features_data.tsv";
output_file	<- "data/target_map/descriptor_all_all_adverse_event_select_features_data_cluster" 
N_repeat	<- 20;

## 1. Read in adverse event ~ target mapping data 
ae_target_map <- read.delim(file = "data/target_map/descriptor_all_all_adverse_event_select_features_data.tsv", header = T, sep = "\t");

## 2. Identify clusters of adverse events using pvcluster
ae_clusters <- NULL;
# repeat the analysis by specified number of runs  
for(i in 1:N_repeat){ 
	set.seed(i)
	# hierarchical clustering by pvclust
	ae_clust <- pvclust(t(ae_target_map), method.dist = "binary", method.hclust = "average", nboot = 100);
	# find significant clusters
	ae_clusters[[i]] <- pvpick(ae_clust)[[1]];	
}

## 3. Build output dataframes of identified clusters 
# collect clusters identified from multiple runs to find unique clusters 
ae_clusters <- unlist(ae_clusters, recursive = F);
ae_clusters <- sapply(ae_clusters, function(ac) paste(sort(ac), collapse = ","));
ae_clusters_table <- table(ae_clusters);
# obtain number of adverse events in each unique cluster 
act_len <- sapply(names(ae_clusters_table), function(nact) length(strsplit(nact, ",")[[1]]));
# build output data frame 
ae_cluster_df <- data.frame(names(ae_clusters_table), act_len, as.numeric(ae_clusters_table));
colnames(ae_cluster_df) <- c("adverse_events_cluster", "N_adverse_events", "N_reproducibility");
# sort rows by number of appearances from multiple runs  
acd_od <- order(ae_clusters_table, decreasing = T);
ae_cluster_df <- ae_cluster_df[acd_od, ];
write.table(ae_cluster_df, paste(output_file, "_ae.tsv", sep = ""), sep = "\t", row.names = F, col.names= T, quote = F);
