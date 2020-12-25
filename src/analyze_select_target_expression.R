# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script analyzes whether selected target genes of adverse events are more likely to be differentially expressed in the matched tissues using GTEx expression data.


## functions
source("src/functions.R");


## 0. Input arguments 
diff_exp_file	<- "data/target_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_uniprot_median_tpm_adjust.tsv";
ae_map_file	<- "data/target_expression/adverse_event_tissue_map.tsv";
select_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_all_adverse_event_feature_select_implementation/descriptor_all_all_adverse_event_select_features.tsv"
out_file	<- "data/target_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_uniprot_median_tpm_adjust_all_adverse_event_select_targets";

## 1. Read in tissue-specific differential expression data, tissue map, and selected targets of adverse events 
# read in the processed GTEx differential expression data
diff_exp <- read.delim(file = diff_exp_file, sep = "\t", header = T); 
# read in adverse event ~ tissue map  
ae_tissue <- read.delim(file = ae_map_file, sep = "\t", header = T);
all_aes <- as.character(ae_tissue$adverse_event);
aa_tissues <- as.character(ae_tissue$tissue_name);
# read in data of selected features
ae_select_features <- read.delim(file = select_file, sep = "\t", header = T); 
# obtain the selected features of each adverse event of interest  
aa_id <- sapply(all_aes, function(aa) which(ae_select_features$adverse_event %in% aa));
aa_features <- lapply(ae_select_features$select_features[aa_id], function(asfsf) strsplit(asfsf, ",", fixed = T)[[1]]);

## 2. Extract differential expression values of selected target genes, compare them with all measured genes 
# iterate by adverse event 
select_diff_exp <- mapply(function(aa, at, af){
	# obtain differential expression values of selected target genes in the specified tissue 
	af_id <- which(rownames(diff_exp) %in% af);
	af_diff_exp <- diff_exp[af_id, at];
	# combine it with differential expression values of all measured genes 
	at_diff_exp <- diff_exp[, at];
	combine_diff_exp <- c(af_diff_exp, at_diff_exp);
	# test whether selected target genes are more likely to be differentially expressed compared to all measure genes using wilcoxon test
	aa_pv <- wilcox.test(af_diff_exp, at_diff_exp, alternative = "greater")$p.value;
	# output the extracted differential expression values along with other info (adverse event term, tissue, gene name)
	combine_genes <- c(rownames(diff_exp)[af_id], rownames(diff_exp));
	combine_type <- c(rep("Selected target genes", length(af_id)), rep("All genes measured by GTEx", nrow(diff_exp)));
	combine_len <- length(combine_diff_exp);
	aa_vec <- rep(aa, combine_len);
	at_vec <- rep(at, combine_len);
	aa_df <- data.frame(aa_vec, at_vec, combine_genes, combine_type, combine_diff_exp);
	colnames(aa_df) <- c("adverse_event", "tissue", "gene", "gene_group", "differential_expression");
	return(ls = list(exp_df = aa_df, exp_pv = aa_pv));	
}, all_aes, aa_tissues, aa_features, SIMPLIFY = F);
# output extracted differential expression values of selected target genes from all adverse events in data frame
select_de <- lapply(select_diff_exp, function(sde) sde[[1]]);
select_de_df <- do.call(rbind, select_de);
write.table(select_de_df, file = paste(out_file, "_compare_de.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# perform multiple testing correction on all adverse events of interest, output results in data frame 
select_de_pv <- sapply(select_diff_exp, function(sde) sde[[2]]);
select_de_fdr <- p.adjust(select_de_pv, method = "fdr"); 
select_de_test_df <- data.frame(all_aes, aa_tissues, select_de_pv, select_de_fdr);
colnames(select_de_test_df) <- c("adverse_event", "tissue", "p_value", "fdr");
write.table(select_de_test_df, file = paste(out_file, "_compare_pv.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);

