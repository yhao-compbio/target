## !/usr/bin/env python
## created by Yun Hao @MooreLab 2019
## This script integrates three different resources to build human druggable genome 

## 1. Obtain entrez-uniprot ID mapping of human genes 
# read in uniprot ID mapping of human proteins 
uniprot_df <- read.delim(file = "downloads/id_map/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+review--.tab", header = T, sep = "\t");
# remove proteins that are not reviewed 
review_id <- which(as.character(uniprot_df$Status) %in% "reviewed");
# obtain set of unique uniprot protein IDs 
review_uniprot <- unique(as.character(uniprot_df$Entry)[review_id]);
# read in entrez ID mapping of human genes  
id_map_df <- read.delim(file = "downloads/id_map/HUMAN_9606_idmapping_selected.tab", header = F, sep = "\t");
id_map_df <- id_map_df[, c(1, 3)];
colnames(id_map_df) <- c("uniprot_id", "GeneID");
# remove genes without proper entrez identifiers 
nna_id <- which(as.character(id_map_df$GeneID) != "");
id_map_df <- id_map_df[nna_id,];
# remove genes with unreviewed protein products 
review_id <- which(as.character(id_map_df$uniprot_id) %in% review_uniprot);
id_map_df <- id_map_df[review_id,];
# separate entrez identifiers when one uniprot protein is mapped to multiple entrez genes
gene_id_list <- lapply(as.character(id_map_df$GeneID), function(imd){
        entrezs <- strsplit(imd, "; ", fixed = T)[[1]];
        return(as.integer(entrezs));
});
gene_id_len <- sapply(gene_id_list, length);
uniprot_id_list <- mapply(function(x, y) rep(x, y), as.character(id_map_df$uniprot_id), gene_id_len);
# obtain full mapping between entrez IDs of genes and uniprot IDs of protein products 
id_map_df1 <- data.frame(unlist(uniprot_id_list), unlist(gene_id_list));
colnames(id_map_df1) <- c("uniprot_id", "GeneID");

## 2. Process protein-class mapping data from dGRene
# read in dGene data
dgene_df <- read.csv("downloads/dgene/journal.pone.0067980.s002.CSV");
dgene_df <- dgene_df[, c("GeneID", "class")];
# map proteins from dGene to their uniprot IDs 
dgene_uniprot_df <- merge(dgene_df, id_map_df1, by = "GeneID");
# standardize protein subclass names 
class_name1 <- c("GPCR", "NHR", "PI3K", "Protease", "PI", "PTEN", "PTP", "MTMR", "STK", "RTK");
names(class_name1) <- c("GPCR", "NHR", "PI3K", "PROTEASE", "PROT_INHIB", "PTEN", "PTP", "PTP_MTMR", "ST_KINASE", "Y_KINASE");
dgene_uniprot_df$class <- class_name1[as.character(dgene_uniprot_df$class)];
dgene_uniprot_df$uniprot_id <- as.character(dgene_uniprot_df$uniprot_id);

## 3. Process protein-class mapping data from GtoPDB 
# read in GtoPDB data 
gtopdb_df <- read.csv(file = "downloads/gtopdb/targets_and_families.csv");
# remove proteins without proper 
nna_id <- which(as.character(gtopdb_df$Human.SwissProt) != "");
gtopdb_df <- gtopdb_df[nna_id,];
# remove proteins that have already been mapped through dGene 
dgene_gtop_id <- which(as.character(gtopdb_df$Human.SwissProt) %in% dgene_uniprot_df$uniprot_id);
other_id <- setdiff(1:nrow(gtopdb_df), dgene_gtop_id);
gtopdb_df <- gtopdb_df[other_id, c("Human.SwissProt", "Type")];
gtopdb_df <- unique(gtopdb_df);
# separate uniprot identifiers when one protein is mapped to multiple uniprot identifiers
uniprot_id_list <- lapply(as.character(gtopdb_df$Human.SwissProt), function(gdhs){
	ids <- strsplit(gdhs, "|", fixed = T)[[1]];
	return(ids);
});
uniprot_id_len <- sapply(uniprot_id_list, length);
type_list <- mapply(function(x,y) rep(x,y), as.character(gtopdb_df$Type), uniprot_id_len);
# obtain full protein-class mapping from GtoPDB 
gtopdb_df1 <- data.frame(unlist(uniprot_id_list), unlist(type_list));
colnames(gtopdb_df1) <- c("uniprot_id", "class");
# map proteins from dGene to their entrez IDs 
gtopdb_entrez_df <- merge(gtopdb_df1, id_map_df1, by = "uniprot_id");
# standardize protein subclass names 
class_name2 <- c("GPCR", "NHR", "VGIC", "LGIC", "OIC", "Transporter", "OCR", "OE", "OP");
names(class_name2) <- c("gpcr", "nhr", "vgic", "lgic", "other_ic", "transporter", "catalytic_receptor", "enzyme", "other_protein");
gtopdb_entrez_df$class <- class_name2[as.character(gtopdb_entrez_df$class)];
gtopdb_entrez_df$uniprot_id <- as.character(gtopdb_entrez_df$uniprot_id);
gtopdb_entrez_df <- gtopdb_entrez_df[,colnames(dgene_uniprot_df)];
# combine mapping results from dGene and GtoPDB
dgene_gtop_df <- rbind(dgene_uniprot_df, gtopdb_entrez_df);

## 4. Process protein-class mapping data from DrugBank 
# Read in drugbank data
drugbank_df <- read.csv(file = "downloads/drugbank/uniprot\ links.csv");
# find proteins that have not been mapped via GtoPDB or dGene 
drugbank_target <- unique(as.character(as.character(drugbank_df$UniProt.ID)));
other_target <- setdiff(drugbank_target, dgene_gtop_df$uniprot_id);
other_target_id <- which(as.character(id_map_df1$uniprot_id) %in% other_target);
# assign these proteins into the group of OP('Other Protein')
other_target_df <- cbind(id_map_df1[other_target_id, ], "OP");
colnames(other_target_df) <- c("uniprot_id", "GeneID", "class");
other_target_df <- other_target_df[, colnames(dgene_gtop_df)];

## 5. Combine mapping results from all three sources
druggable_df <- rbind(dgene_gtop_df, other_target_df);
# group subclasses into general functional classes 
cate_name <- c("G protein-coupled receptor", "Nuclear hormone receptor", "Ion channel", "Ion channel", "Ion channel", "Transporter", "Catalytic receptor", "Catalytic receptor", "Catalytic receptor", "Enzyme", "Enzyme", "Enzyme", "Enzyme", "Enzyme", "Enzyme", "Enzyme", "Other protein");
names(cate_name) <- c("GPCR", "NHR", "VGIC", "LGIC", "OIC", "Transporter", "STK", "RTK", "OCR", "PI3K", "Protease", "PI", "PTEN", "PTP", "MTMR", "OE", "OP");
druggable_df$cate <- cate_name[druggable_df$class];
colnames(druggable_df) <- c("entrez_id", "subclass", "uniprot_id", "class");
druggable_df <- druggable_df[,c("uniprot_id", "entrez_id", "class", "subclass")];
# sort proteins by class names, then subclass names 
od <- order(druggable_df$class, druggable_df$subclass);
druggable_df <- druggable_df[od, ];
# output the final mapping results  
write.table(druggable_df, file = "data/human_druggable_genome_protein_class.tsv", sep = "\t", col.names = T, row.names = F, quote = F);
