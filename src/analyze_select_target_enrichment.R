# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This scripts performs Gene Ontology/Reactome pathway enrichment analysis of target genes of interest.


## functions
library(GO.db);
source("src/functions.R");


## 0. Input arguments
select_feature_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_all_adverse_event_feature_select_implementation/descriptor_all_all_adverse_event_select_features.tsv";
ae_file			<- "data/target_enrichment/adverse_event_of_interest.txt";
reactome_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";
go_file			<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/gene_ontology/goa_human.gaf";
min_path_size		<- 10;
max_path_size		<- 100;
output_file		<- "data/target_enrichment/descriptor_all_all_adverse_event_select_features";

## 1. Obtain the selected features of adverse events 
# read in data frame of selected features
ae_select_feature_df <- read.delim(file = select_feature_file, sep = "\t", header = T);
# read in adverse event of interest 
all_aes <- readLines(ae_file);
aa_id <- sapply(all_aes, function(aa) which(ae_select_feature_df$adverse_event %in% aa));
ae_select_feature_df <- ae_select_feature_df[aa_id, ];
# obtain the selected features of each adverse event
ae_select_features <- lapply(ae_select_feature_df$select_features, function(asfsf) strsplit(asfsf, ",", fixed = T)[[1]]);
names(ae_select_features) <- sapply(ae_select_feature_df$adverse_event, function(asfae) paste(strsplit(asfae, "_")[[1]], collapse = " ")); 

## 2. Perform Reactome pathway enrichment analysis on selected target genes   
# read in reactome annotation data 
reactome_df <- read.delim(file = reactome_file, header = F, sep = "\t");
# remove non-human annotations
human_id <- which(reactome_df[,6] %in% "Homo sapiens");
reactome_df <- reactome_df[human_id, c(1, 4)];
# perform Reactome pathway enrichment
ae_reactome_enrich <- mapply(function(asf, nasf){
	asf_enrich <- performm.enrichment.analysis(asf, reactome_df, min_path_size, max_path_size, "fdr");
	adverse_event <- rep(nasf, nrow(asf_enrich));
	asf_df <- cbind(adverse_event, asf_enrich);
        return(asf_df);
}, ae_select_features, names(ae_select_features), SIMPLIFY = F);
# combine enrichment analysis results of all adverse events, output in data frame  
ae_reactome_enrich_df <- do.call(rbind, ae_reactome_enrich);
write.table(ae_reactome_enrich_df, file = paste(output_file, "_reactome_enrichment_results.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);

## 3. Perform Reactome pathway enrichment analysis on selected target genes    
# read in human GO annotation data
go_df <- read.delim(file = go_file, header = F, sep = "\t", skip = 38);
# map GO term IDs to their function descriptions  
go_info <- as.list(GOTERM);
go_id <- sapply(go_info, function(gi) GOID(gi));
go_term <- sapply(go_info, function(gi) Term(gi));
go_id_term <- data.frame(go_id, go_term);
go_term_df <- merge(go_df, go_id_term, by.x = "V5", by.y = "go_id");
# group GO annotations into BP (biological process), MF (molecular function), and CC (cellular component)
bp_id <- which(go_term_df$V9 == "P");
bp_term_df <- go_term_df[bp_id, c(3, 18)];
mf_id <- which(go_term_df$V9 == "F");
mf_term_df <- go_term_df[mf_id, c(3, 18)];
cc_id <- which(go_term_df$V9 == "C");
cc_term_df <- go_term_df[cc_id, c(3, 18)];
# perform GO BP term enrichment 
ae_bp_enrich <- mapply(function(asf, nasf){
	asf_enrich <- performm.enrichment.analysis(asf, bp_term_df, min_path_size, max_path_size, "fdr");
	adverse_event <- rep(nasf, nrow(asf_enrich));
	asf_df <- cbind(adverse_event, asf_enrich);
	return(asf_df);
}, ae_select_features, names(ae_select_features), SIMPLIFY = F);
# combine enrichment analysis results of all adverse events, output in data frame 
ae_bp_enrich_df <- do.call(rbind, ae_bp_enrich);
write.table(ae_bp_enrich_df, file = paste(output_file, "_go_bp_enrichment_results.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# perform GO MF term enrichment 
ae_mf_enrich <- mapply(function(asf, nasf){
	asf_enrich <- performm.enrichment.analysis(asf, mf_term_df, 10, 100, "fdr");
	adverse_event <- rep(nasf, nrow(asf_enrich));
	asf_df <- cbind(adverse_event, asf_enrich);
	return(asf_df);
}, ae_select_features, names(ae_select_features), SIMPLIFY = F);
# combine enrichment analysis results of all adverse events, output in data frame  
ae_mf_enrich_df <- do.call(rbind, ae_mf_enrich);
write.table(ae_mf_enrich_df, file = paste(output_file, "_go_mf_enrichment_results.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# perform GO CC term enrichment 
ae_cc_enrich <- mapply(function(asf, nasf){
	asf_enrich <- performm.enrichment.analysis(asf, cc_term_df, 10, 100, "fdr");
	adverse_event <- rep(nasf, nrow(asf_enrich));
	asf_df <- cbind(adverse_event, asf_enrich);
	return(asf_df);
}, ae_select_features, names(ae_select_features), SIMPLIFY = F);
# combine enrichment analysis results of all adverse events, output in data frame 
ae_cc_enrich_df <- do.call(rbind, ae_cc_enrich);
write.table(ae_cc_enrich_df, file = paste(output_file, "_go_cc_enrichment_results.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
