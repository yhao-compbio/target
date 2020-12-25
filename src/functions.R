# !/usr/bin/env Rscript
# created by Yun Hao @MooreLab 2019
# This script contains R functions required for other scripts in the repository.


## This function groups a vector by categories of its elements.
group.vector.by.categories <- function(cate, vec){
	## 0. Input arguments 
		# cate: category of vectors  
		# vec: vector
	
	## 1. Sort vector by order of categories
	vec_od <- order(cate);
	cate <- cate[vec_od];
	vec <- vec[vec_od];

	## 2. Group elements of same category together
	# obtain unique categories
	cate_table <- table(cate);
	# obtain lower bound index/upper bound index of each unique category
	lower_ids <- cumsum(cate_table);
	upper_ids <- c(0, lower_ids[-length(lower_ids)]) + 1;
	# return list of vectors 
	vec_list <- mapply(function(li, ui) vec[li:ui], lower_ids, upper_ids, SIMPLIFY=F);
	names(vec_list) <- names(cate_table);

	return(vec_list);
}


## This function merges two lists into one. 
merge.list <- function(l1, l2){
	## 0. Input arguments: 
		# l1: first list to be merged 
		# l2: second list to be merged

	## 1. Obtain all unique list names  
	l_union <- unique(union(names(l1), names(l2)));
	
	## 2. Obtain elements of each list name 
	# iterate by list name 
	m_list <- lapply(l_union,function(x){
		# obtain name index in the two lists, respectively   
		l1_id <- which(names(l1) %in% x)
		l2_id <- which(names(l2) %in% x)
		# combine elements if the name appears in both list  
		if(length(l1_id) > 0 && length(l2_id) > 0){
			element <- c(unlist(l1[l1_id]),unlist(l2[l2_id]))
		}
		# get elements from the list if the name only appaers in one list
                else if(length(l1_id) == 0 && length(l2_id) > 0){
			element <- c(unlist(l2[l2_id]))
		} else{
			element <- c(unlist(l1[l1_id]))
		}
		element <- unique(element);
		return(element);
	});
	names(m_list) <- l_union;
        return(m_list);
}


## This function performs Fisher Exact Test between a set of interesting genes and a list of gene sets.
compute.fisher.p.value <- function(geneSetList, interestingGene, background){
	## 0. Input arguments: 
		# geneSetList: list of gene sets
	 	# interestingGene: interesting genes 
		# background: background genes

	## 1. Remove genes that are not in the background list from the ones to be tested    
	interestingGene <- intersect(interestingGene,background);
	geneSetList <- lapply(geneSetList,function(x) intersect(x,background));

	## 2. Perform Fisher Exact Test for each gene set  
	# iterate by gene set 
	fisher_pval <- sapply(geneSetList,function(gs){
		# count numbers in the 2*2 contigency table  
		GsIn <- length(intersect(gs,interestingGene));
		diffin <- setdiff(background,interestingGene);
		diffgs <- setdiff(background,gs);
		DiffgsIn <- length(intersect(interestingGene,diffgs));
		DiffinGs <- length(intersect(gs,diffin));
		DiffinDiffgs <- length(intersect(diffin,diffgs));
		# compute P value of Fisher Exact Test  
		fisher.test(matrix(c(GsIn,DiffgsIn,DiffinGs,DiffinDiffgs),2,2),alternative="greater")$p.value;
	});
	names(fisher_pval) <- names(geneSetList);
	return(fisher_pval);
}


## This function performs gene set enrichment analysis between a set of interesting genes and a list of gene sets. 
performm.enrichment.analysis <- function(query_genes, annotation_df, min_size, max_size, adjust_method, only_sig = TRUE, sig_threshold = 0.05){
        ## 0. Input arguments
                # query_genes: interesting genes
                # annotation_df: data frame that contains gene-gene set annotation 
		# min_size: minimal size for a gene set to be considered 
                # max_size: maximal size for a gene set to be considered 
		# adjust_method: method to perform multiple testing correction
		# only_sig: whether to only return the results of significant pathways (default: TRUE) 
		# sig_threshold: adjusted p value threshold of significance 

        ## 1. Obtain gene sets from annotation data 
	# group gene annotations by gene sets  
	path_genes <- group.vector.by.categories(annotation_df[,2], annotation_df[,1]);
	# obtain all genes with annotations as background genes 
	bg_genes <- unique(annotation_df[,1]);
	# remove gene sets that are not in the specified range
	path_len <- sapply(path_genes, length);
	path_id <- intersect(which(path_len >= min_size), which(path_len <= max_size));
	path_genes <- path_genes[path_id];

	## 2. Perform enrichment analysis 
	# compute p value for each pathway  
	path_pv <- compute.fisher.p.value(path_genes, query_genes, bg_genes);
	path_pv <- sort(path_pv);
	# adjust p value by the specified method  
	path_adj_p <- p.adjust(path_pv, method = adjust_method);
	# remove results of insignificant pathways, if specified so
	if(only_sig == TRUE){
		sig_id <- which(path_adj_p < sig_threshold)
		path_pv <- path_pv[sig_id]
		path_adj_p <- path_adj_p[sig_id]
	}
	path_pv_len <- path_len[names(path_pv)];

	## 3. Build data frame that shows the result of enrichment analysis 
	result_df <- data.frame(names(path_pv), path_pv_len, path_pv, path_adj_p);
	colnames(result_df) <- c("pathway", "size", "p_value", "adjusted_p_value");
	return(result_df);
}


## This function computes the similarity of functional annotations between input target pairs.
compute.targets.function.similarity <- function(target1, target2, target_annotation_df, min_size, max_size){
	## 0. Input arguments
		# target1: vector that contains the first target of the pair 
		# target2: vector that contains the second target of the pair 
		# target_annotation_df: data frame that contains the function annotation of targets 
		# min_size: minimal size for a function set to be considered
		# max_size: maximal size for a function set to be considered

	## 1. Obtain the functional annotation of each target  
	# remove function sets that are not in the specified range
	path_len <- table(as.character(target_annotation_df[,2]));
	len_id <- intersect(which(path_len >= min_size), which(path_len <= max_size));
	all_paths <- names(path_len)[len_id];
	# remove the function annotations of the above sets 
	path_id <- which(as.character(target_annotation_df[,2]) %in% all_paths);
	target_annotation_df <- target_annotation_df[path_id, ];
	# remove the function annotations of targets that are not in the pairs  
	all_targets <- unique(c(target1, target2));	
	target_id <- which(as.character(target_annotation_df[,1]) %in% all_targets);
	target_annotation_df <- target_annotation_df[target_id, ];
	# group annotation sets by targets  
	target_annotation <- group.vector.by.categories(as.character(target_annotation_df[,1]), as.character(target_annotation_df[,2]));
	
	## 2. Compute the similarity of functional annotations between target pairs 
	target_sim <- mapply(function(t1, t2){
		# obtain target index in the annotation list  
		t1_id <- which(names(target_annotation) %in% t1);
		t2_id <- which(names(target_annotation) %in% t2);
		# if both targets have annotations, compute jaccard similarity of annotation sets 
		if ((length(t1_id) > 0) && (length(t2_id) > 0)){
			N_inter <- length(intersect(target_annotation[[t1_id]], target_annotation[[t2_id]]))
			N_union <- length(union(target_annotation[[t1_id]], target_annotation[[t2_id]]))
			sim <- N_inter/N_union
		}
		# otherwise, return NA as results   
		else{
			sim <- NA
		}
		return(sim);			
	}, target1, target2);
	sim_df <- data.frame(target1, target2, round(target_sim, 5));
	colnames(sim_df) <- c("target1" , "target2", "function_similarity");

	return(sim_df);
}


## This function groups target pairs into the ones with common funtions and the ones without common functions, the return along with their similarity values of selected features.  
generate.function.group.dataframe <- function(func_sim_df, feat_sim_df, annot_source){
	## 0. Input arguments: 
		# func_sim_df: data frame that contains function similarity of target pairs  
		# feat_sim_df: data frame that contains function similarity of target pairs 
		# annot_source: name of the annotation source by which the function similarity was computed 
	 
	## 1. Build data frame that contains the feature similarity of target pairs with common funtions 
	# obtain target pairs with common funtions (function similarity > 0)
 	common_id <- which(func_sim_df$function_similarity > 0);
	N_common <- length(common_id);
	# obtain the feature similarity of these pairs 
	func_group <- rep("pairs w common function", N_common);
	common_sim_df <- data.frame(func_sim_df$target1[common_id], func_sim_df$target2[common_id], feat_sim_df$feature_jaccard_similarity[common_id], func_group);
	colnames(common_sim_df) <- c("target1", "target2", "feature_similarity", "function_similarity_group");
	
	## 2. Build data frame that contains the feature similarity of target pairs without common funtions 
	# obtain target pairs without common funtions (function similarity = 0)
	no_common_id <- which(func_sim_df$function_similarity == 0);
	N_no_common <- length(no_common_id);
	# obtain the feature similarity of these pairs 
	func_group <- rep("pairs w/o common function", N_no_common);
	no_common_sim_df <- data.frame(func_sim_df$target1[no_common_id], func_sim_df$target2[no_common_id], feat_sim_df$feature_jaccard_similarity[no_common_id], func_group);
	colnames(no_common_sim_df) <- c("target1", "target2", "feature_similarity", "function_similarity_group");
		
	## 3. Combine the two data frames above
	func_group_sim_df <- rbind(common_sim_df, no_common_sim_df);
	func_group_sim_df$annotation_source <- rep(annot_source, nrow(func_group_sim_df));

	return(func_group_sim_df);
}
