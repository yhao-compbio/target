# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script analyzes whether selected targets are more likely to be disease-related genes using CTD gene-diseaese connection data.


## functions
source("src/functions.R");


## 0. Input arguments
id_map_file	<- "downloads/id_map/hgnc_gene_names.tsv";
bg_file		<- "data/target_disease/CTD_genes_diseases_direct_evidence_genes.txt";
connection_file	<- "data/target_disease/CTD_genes_diseases_direct_evidence_diseases_filtered_connection_map.tsv"
select_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_all_adverse_event_feature_select_implementation/descriptor_all_all_adverse_event_select_features.tsv";
map_file	<- "data/target_disease/adverse_event_disease_term_map.tsv";
output_file	<- "data/target_disease/descriptor_all_all_adverse_event_select_features";

## 1. Obtain full map bewteen Uniprot IDs and symbol IDs 
# read in uniprot-to-symbol mapping data from HGNC 
id_map <- read.delim(file = id_map_file, header = T, sep = "\t");
# remove gene symbols that cannot be mapped to a UniProt ID 
matched_id <- which(id_map$UniProt.ID.supplied.by.UniProt. != "");
id_map <- id_map[matched_id, ];
# obtain mapped symbol IDs of each Uniprot ID 
id_uni_ids <- lapply(id_map$UniProt.ID.supplied.by.UniProt., function(idu) strsplit(idu, ", ")[[1]]);
id_uni_len <- sapply(id_uni_ids, length);
# build full map of Uniprot ID ~ symbol ID
id_uni_sym <- mapply(function(imas, iul){
	rep(imas, iul)
}, id_map$Approved.symbol, id_uni_len, SIMPLIFY = F);
id_full_map_df <- data.frame(unlist(id_uni_sym), unlist(id_uni_ids));
colnames(id_full_map_df) <- c("Approved.symbol", "UniProt.ID.supplied.by.UniProt.");

## 2. Map selected target genes to their symbol IDs
# read in data frame of selected features 
ae_select_features_df <- read.delim(file = select_file, sep = "\t", header = T);
# obtain the selected features of each adverse event of interest  
ae_select_features <- lapply(ae_select_features_df$select_features, function(asfdsf) strsplit(asfdsf, ",")[[1]]);
names(ae_select_features) <- ae_select_features_df$adverse_event;
# map Uniprot IDs of selected features to their symbol IDs
ae_select_targets <- mapply(function(asf, nasf, asfmddg){
	# iterate by adverse event
	asf_sym <- sapply(asf, function(a){
		a_id <- which(id_full_map_df$UniProt.ID.supplied.by.UniProt. %in% a);
		if(length(a_id) == 0)	return(NA)
		else	return(id_full_map_df$Approved.symbol[[a_id]])
	});
	# combine mapped symbol IDs along with other info (adverse event, disease group, Uniprot ID) 
	asf_df <- data.frame(rep(nasf, length(asf)), asf, asf_sym);
	colnames(asf_df) <- c("adverse_event", "select_target_uniprot", "select_target_symbol");	
	return(asf_df);
}, ae_select_features, names(ae_select_features), SIMPLIFY= F);
# combine mapped results from all adverse events, output in data frame 
all_ae_select_targets_df <- do.call(rbind, ae_select_targets);
write.table(all_ae_select_targets_df, file = paste(output_file, "_symbol_all.tsv", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F);

## 3. 
# read in adverse event ~ disease group map 
ae_map_df <- read.delim(file = map_file, sep = "\t", header = T);
# map selected targets of adverse events to their associated disease groups 
ae_select_targets_df <- merge(all_ae_select_targets_df, ae_map_df, by = "adverse_event");
ae_select_targets_df <- ae_select_targets_df[,c("adverse_event", "disease_group", "select_target_uniprot", "select_target_symbol")];
ae_select_targets_od <- order(ae_select_targets_df$disease_group, ae_select_targets_df$adverse_event, ae_select_targets_df$select_target_uniprot);
ae_select_targets_df <- ae_select_targets_df[ae_select_targets_od, ];
write.table(all_ae_select_targets_df, file = paste(output_file, "_symbol.tsv", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F);

## 3. Connect selected targets to associated diseases with direct evidence supporting  
# read in CTD gene-disease connection data    
connection_df <- read.delim(file = connection_file, header = T, sep = "\t");
# identify associated diseases of selected targets, output in data frame along with other info   
ae_select_targets_df$select_target_disease <- mapply(function(astddg, astdsts){
	# iterate by target gene-adverse event, obtain the subset of disease terms that belong to the diseaese group of the adverse event 
	astddg_id <- which(connection_df$disease_group %in% astddg);
	astddg_df <- connection_df[astddg_id, ];
	# identify the disease terms connected to the target gene  
	astdsts_id <- which(astddg_df$target_symbol %in% astdsts);
	if(length(astdsts_id) == 0)	return(NA)
	else{
		astdsts_terms <- astddg_df$disease_term[astdsts_id]
		astdsts_term_char <- paste(astdsts_terms, collapse = ";")
		return(astdsts_term_char)
	}
}, ae_select_targets_df$disease_group, ae_select_targets_df$select_target_symbol);
write.table(ae_select_targets_df, file = paste(output_file, "_symbol_connection_map.tsv", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F);

# 4. Perform disease gene enrichment analysis on selected targets of adverse events  
# read in list of all genes that are connected to at least one disease with direct evidence (background genes)  
all_direct_genes <- readLines(bg_file);
# perform disease gene enrichment analysis on each set of selected targets 
target_disease_enrich <- mapply(function(asfmdae, asfmddg){
	# iterate by adverse event, obtain the selected targets for the adverse event 
	asfmdae_id <- which(ae_select_targets_df$adverse_event %in% asfmdae);
	asfmdae_df <- ae_select_targets_df[asfmdae_id, ];
	# remove selected targets that cannot be mapped to a symbol ID
	nna_id1 <- which(is.na(asfmdae_df$select_target_symbol) == F);
	asfmdae_df <- asfmdae_df[nna_id1, ];
	# identify subset of selected target genes connected to associated disease terms by direct evidence, and its complement set in the background genes
	asfmdae_targets <- asfmdae_df$select_target_symbol;
	asfmdae_targets <- intersect(asfmdae_targets, all_direct_genes);
	asfmdae_targets_c <- setdiff(all_direct_genes, asfmdae_targets);
	# identify all genes connected to associated disease terms by direct evidence, and its complement set in the background genes
	asfmddg_id <- which(connection_df$disease_group %in% asfmddg);
	asfmddg_targets <- unique(connection_df$target_symbol[asfmddg_id]);
	asfmddg_targets_c <- setdiff(all_direct_genes, asfmddg_targets);
	# count the numbers in 2*2 contigency table  
	m11 <- length(which(is.na(asfmdae_df$select_target_disease) == F));
	m12 <- length(asfmdae_targets) - m11;
	m21 <- length(asfmddg_targets) - m11;
	m22 <- length(intersect(asfmdae_targets_c, asfmddg_targets_c));
	ct22 <- c(m11, m12, m21, m22);
	# compute Fisher Exact test p value 
	ct22_test <- fisher.test(matrix(ct22, 2, 2, byrow = T), alternative = "greater")
	result_vec <- c(ct22, ct22_test$estimate, ct22_test$p.value);
	return(result_vec);
}, ae_map_df$adverse_event, ae_map_df$disease_group);
# combine enrichment analysis results from all adverse events of interest 
target_disease_enrich <- data.frame(t(target_disease_enrich));
colnames(target_disease_enrich) <- c("n_selected_targets_w_direct_evidence", "n_selected_targets_wo_direct_evidence", "n_other_targets_w_direct_evidence", "n_other_targets_wo_direct_evidence", "odds_ratio", "p_value");
# perform multiple testing correction across all adverse events of interest  
target_disease_enrich$fdr <- p.adjust(target_disease_enrich$p_value, method = "fdr");
target_disease_enrich <- round(target_disease_enrich, 3);
# output enrichment analysis result in data frame 
target_disease_enrich1 <- data.frame(ae_map_df$disease_group, ae_map_df$adverse_event, target_disease_enrich);
colnames(target_disease_enrich1) <- c("organ", "adverse_event_term", colnames(target_disease_enrich));
tde_od <- order(target_disease_enrich1$fdr);
target_disease_enrich1 <- target_disease_enrich1[tde_od,];
write.table(target_disease_enrich1, file = paste(output_file, "_symbol_connection_enrich.tsv", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F);

