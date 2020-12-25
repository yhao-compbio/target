# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script processes raw CTD gene-disease connection data for identifying disease terms of interest.  


## functions
source("src/functions.R");


## 0. Input arguments
ctd_file		<- "downloads/ctd/CTD_genes_diseases.tsv.gz";
disease_map_file	<- "data/target_disease/adverse_event_disease_term_map.tsv"
output_file		<- "data/target_disease/CTD_genes_diseases"

## 1. Obtain CTD gene-disease connections with direct evidence 
# read in gene-disease connection data from CTD 
ctd_df <- read.delim(file = "downloads/ctd/CTD_genes_diseases.tsv.gz", header = F, skip = 29);
# identify connections with direct evidence, output these connections 
direct_id <- which(ctd_df$V5 != "");
ctd_df <- ctd_df[direct_id, ];
write.table(ctd_df, file = paste(output_file, "_direct_evidence_data.tsv", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F);
# obtaini all genes that are connected with at least one disease by direct evidence   
direct_ctd_genes <- unique(ctd_df$V1);
writeLines(direct_ctd_genes, paste(output_file, "_direct_evidence_genes.txt", sep = ""));

## 2. Obtain gene-disease connections that are associated with disease groups of interest
# read in map of disease group ~ search key words  
disease_map_df <- read.delim(file = disease_map_file, header = T, sep = "\t");
disease_map_df <- unique(disease_map_df[, c("disease_group", "disease_term")]);
# obtain search key words for each disease group 
group_terms <- lapply(disease_map_df$disease_term, function(dmddt) strsplit(dmddt, ",", fixed = T)[[1]]);
names(group_terms) <- disease_map_df$disease_group;
# implement string match between search key words and all disease terms 
ctd_df$V3 <- tolower(ctd_df$V3);
all_diseases <- unique(ctd_df$V3);
group_diseases <- mapply(function(gt, ngt){
	# identify disease terms that contain the key words 
	gt_id <- lapply(gt, function(tt){
		tt_id <- grep(tt, all_diseases);
		return(tt_id);
	});
	all_gt_ids <- unique(unlist(gt_id));
	# output disease group and terms in data frame   
	gt_df <- data.frame(rep(ngt, length(all_gt_ids)), all_diseases[all_gt_ids]);
	colnames(gt_df) <- c("group", "disease_term");
	return((gt_df));
}, group_terms, names(group_terms), SIMPLIFY = F);
# output identified disease terms of all groups 
group_disease_df <- do.call(rbind, group_diseases);
write.table(group_disease_df, file = paste(output_file, "_direct_evidence_diseases.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);

## 3. Manually review the matched disease terms of disease groups. Remove terms that have no associations with the disease group (Some unrelated terms can be identified via string match (i.e. 'language delay and attention deficit-hyperactivity disorder/cognitive impairment with or without cardiac arrhythmia' matched via search of 'cardiac'). 
# manually review, store terms to be moved in ""_direct_evidence_diseases_remove.tsv"".
# read in terms to be removed   
disease_remove_df <- read.delim(file = paste(output_file, "_direct_evidence_diseases_remove.tsv", sep = ""), header = T, sep = "\t");
# group terms to be removed by disease group 
disease_remove <- group.vector.by.categories(disease_remove_df$disease_group, disease_remove_df$disease_term);
# remove the terms by group  
group_diseases1 <- mapply(function(gd, ngd){
	ngd_id <- which(names(disease_remove) %in% ngd);
	if(length(ngd_id) == 0)	return(gd)
	else{
		gd_terms <- setdiff(gd$disease_term, disease_remove[[ngd_id]])
		new_gd <- data.frame(rep(ngd, length(gd_terms)), gd_terms)
		colnames(new_gd) <- c("group", "disease_term")
		return(new_gd)
	}		
}, group_diseases, names(group_diseases), SIMPLIFY = F);
# output the filtered terms 
group_disease_df1 <- do.call(rbind, group_diseases1);
write.table(group_disease_df1, file = paste(output_file, "_direct_evidence_diseases_filtered.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);

## 4. Map filtered disease terms to connected genes with direct evidence, output in data frame along with other info  
ctd_disease_df <- merge(ctd_df, group_disease_df1, by.x = "V3", by.y = "disease_term");
ctd_disease_out_df <- ctd_disease_df[, c(2, 10, 1, 5)];
ctd_disease_out_df <- unique(ctd_disease_out_df);
colnames(ctd_disease_out_df) <- c("target_symbol", "disease_group", "disease_term", "connection_evidence");
write.table(ctd_disease_out_df, file = paste(output_file, "_direct_evidence_diseases_filtered_connection_map.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
