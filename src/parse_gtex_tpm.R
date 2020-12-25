# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2020
# This script processes raw GTEx tissue-specific expression data for gene ID mapping and differential expression analysis. 


## functions
source("src/functions.R");


## 0. Input arguments
gtex_file	<- "downloads/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz";
id_map_file	<- "downloads/id_map/hgnc_gene_names.tsv";
output_file	<- "data/target_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9";

## 1. Obtain GTEx tissue-specific expression data 
# read in GTEx tissue-specific expression data 
gtex_tpm <- read.delim(file = gtex_file, header = T, sep = "\t", skip = 2);
# obtain GTEx tissue names 
tissue_id <- setdiff(1:ncol(gtex_tpm), which(colnames(gtex_tpm) %in% c("Name", "Description")));
tissue_name <- sapply(colnames(gtex_tpm)[tissue_id], function(cgt){
	cgt_s <- strsplit(cgt, ".", fixed = T)[[1]];
	cgt_s <- cgt_s[cgt_s != ""];
	return(paste(cgt_s, collapse = "_"));
});
# remove genes that have no expression in all tissues 
tissue_sum <- rowSums(gtex_tpm[ ,tissue_id]);
exp_id <- which(tissue_sum > 0);
gtex_tpm <- gtex_tpm[exp_id, ];

## 2. Obtain mapping bewteen GTEx-measured genes and UniProt IDs 
# read in uniprot-to-symbol mapping data from HGNC 
id_map <- read.delim(file = id_map_file, header = T, sep = "\t");
# identify gene symbols that cannot be mapped to a UniProt ID 
unmatched_id <- which(id_map$UniProt.ID.supplied.by.UniProt. == "");
unmatched_symbol <- unique(id_map$Approved.symbol[unmatched_id]); 
# obtain row index of these genes in the expression matrix 
um_sym_gtex_id <- lapply(unmatched_symbol, function(us){
	us_id <- which(gtex_tpm$Description %in% us);
	return(us_id);
});
names(um_sym_gtex_id) <- unmatched_symbol;
# identify gene symbols that can be mapped to UniProt IDs, group them by unique UniProt ID
matched_id <- which(id_map$UniProt.ID.supplied.by.UniProt. != "");
id_map <- id_map[matched_id, ];
# obtain mapped symbol IDs of each Uniprot ID
id_uni_ids <- lapply(id_map$UniProt.ID.supplied.by.UniProt., function(idu) strsplit(idu, ", ")[[1]]);
id_uni_len <- sapply(id_uni_ids, length);
# group gene symbol IDs by Uniprot IDs
id_uni_sym <- mapply(function(imas, iul){
	rep(imas, iul)
}, id_map$Approved.symbol, id_uni_len, SIMPLIFY = F);
id_full_map_df <- data.frame(unlist(id_uni_sym), unlist(id_uni_ids));
colnames(id_full_map_df) <- c("Approved.symbol", "UniProt.ID.supplied.by.UniProt.");
uni_sym <- group.vector.by.categories(id_full_map_df$UniProt.ID.supplied.by.UniProt., id_full_map_df$Approved.symbol);
# obtain row index of these genes in the expression matrix 
uni_sym_gtex_id <- lapply(uni_sym, function(us){
	us_id <- which(gtex_tpm$Description %in% us);
	return(us_id);
});
# combine the index mapping results of the two groups of genes above
sym_gtex_id <- merge.list(um_sym_gtex_id, uni_sym_gtex_id)
sym_len <- sapply(sym_gtex_id, length);
sym_gtex_id <- sym_gtex_id[sym_len > 0];

## 3. Perform differential expression of genes across all tissues 
# obtain absolute expression values of all genes  
uni_gtex_tpm <- mapply(function(sgi){
	if(length(sgi) == 1)	sgi_exp <- as.numeric(gtex_tpm[sgi, tissue_id])	# when a single row is mapped to the gene ID, take the value of the row
	else	sgi_exp <- as.numeric(colMeans(gtex_tpm[sgi, tissue_id]))	# when multiple rows are mapped to the same ID, take the average all mapped rows
	names(sgi_exp) <- tissue_name;
	return(sgi_exp);
}, sym_gtex_id);
uni_gtex_tpm <- t(uni_gtex_tpm);
# obtain the minimal non-zero expression values in the matrix 
uni_gtex_tpm_vec <- c(as.matrix(uni_gtex_tpm));
tpm_min <- min(uni_gtex_tpm_vec[which(uni_gtex_tpm_vec != 0)]);
# fill zero expression values with 1/10 of the minimal non-zero value (so that the median won't be 0, and be used as denominator in the next step)
uni_gtex_tpm1 <- mapply(function(cugt){
	# find tissues with 0 expression values  
	cugt_exp <- uni_gtex_tpm[ ,cugt]
	z_id <- which(cugt_exp == 0);
	nz_id <- which(cugt_exp != 0);
	# fill 0 with non-negative min/10
	cugt_exp[z_id] <- tpm_min/10;
	return(cugt_exp);
}, colnames(uni_gtex_tpm));
# adjust expression of each gene by the baseline expression (differential expression value)
uni_gtex_tpm1_adj <- mapply(function(rugt){
	# compute median expression of the gene across all tissues 
	rugt_med <- median(uni_gtex_tpm1[rugt, ])
	# divide absolute expression value by the median  
	rugt_tpm_adj <- uni_gtex_tpm1[rugt, ]/rugt_med;
	# take absolute value of log10
	rugt_tpm_adj <- abs(log(rugt_tpm_adj, base = 10));
	return(rugt_tpm_adj);
}, rownames(uni_gtex_tpm1));
uni_gtex_tpm1_adj <- t(uni_gtex_tpm1_adj);
uni_gtex_tpm1_adj <- round(uni_gtex_tpm1_adj, 6);
# output differential expression values 
write.table(uni_gtex_tpm1_adj, file = paste(output_file, "_uniprot_median_tpm_adjust.tsv", sep = ""), sep = "\t", col.names = T, row.names = T, quote = F);
