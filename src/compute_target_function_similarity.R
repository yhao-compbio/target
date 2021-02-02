# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script computes the similarity of functional annotations between given target pairs. 


## functions
source("src/functions.R");


pair_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_0.25_binary_feature_select_implementation/descriptor_all_compare/descriptor_all_select_features_mc_0.85_feature_similarity.tsv"	# name of file that contains the target pairs
output_file	<- "data/function_similarity/descriptor_all_select_features_mc_0.85";     # output file 
reactome_file	<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";
go_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/gene_ontology/goa_human.gaf";

## 1. Obtain query target pairs
# read in input data frame that contains target pairs 
pair_df <- read.delim(file = pair_file, header = T, sep = "\t");
# obtain names of 1st and 2nd targets in the pairs 
pair_target1 <- sapply(as.character(pair_df[,1]), function(acpd) strsplit(acpd, "_", fixed = T)[[1]][[1]]);
pair_target2 <- sapply(as.character(pair_df[,2]), function(acpd) strsplit(acpd, "_", fixed = T)[[1]][[1]]);
all_targets <- unique(c(pair_target1, pair_target2));

## 2. Compute the similarity of Reactome annotations between target pairs  
# read in reactome annotation data 
reactome_annotation_df <- read.delim(file = reactome_file, header = F, sep = "\t");
# remove non-human annotations
human_id <- which(reactome_annotation_df[,6] %in% "Homo sapiens"); 
reactome_annotation_df <- reactome_annotation_df[human_id, c(1,2)];
# compute and output the similarity of Reactome annotations between target pairs   
pair_target_func_sim <- compute.targets.function.similarity(pair_target1, pair_target2, reactome_annotation_df, 10, 100);
write.table(pair_target_func_sim, file = paste(output_file, "_function_similarity_reactome.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# group target pairs into ones with common Reactome funtions and ones withot common Reactome functions
group_sim_reactome_df <- generate.function.group.dataframe(pair_target_func_sim, pair_df, "Reactome");

## 3. Compute the similarity of GO annotations between target pairs  
# read in human GO annotation data 
go_annotation_df <- read.delim(file = go_file, header = F, sep = "\t", skip = 38);
go_annotation_df <- go_annotation_df[, c(2,5)];
# compute and output the similarity of GO annotations between target pairs
pair_target_func_sim <- compute.targets.function.similarity(pair_target1, pair_target2, go_annotation_df, 10, 100);
write.table(pair_target_func_sim, file = paste(output_file, "_function_similarity_go.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# group target pairs into ones with common GO funtions and ones withot common GO functions
group_sim_go_df <- generate.function.group.dataframe(pair_target_func_sim, pair_df, "Gene Ontology");
# combine the two data frame generated above, output along with other info from the input file  
group_sim_df <- rbind(group_sim_reactome_df, group_sim_go_df);
write.table(group_sim_df, file = paste(output_file, "_function_similarity_group.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);

