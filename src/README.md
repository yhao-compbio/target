# This folder contains source code used by the repository.

## R/python scripts 

+ [`integrate_druggable_genome.R`](integrate_druggable_genome.R) integrates three different resources to build human druggable genome. 

+ [`compute_target_similarity.R`](compute_target_similarity.R) computes the similarity of functional annotations between given target pairs.

+ [`visualize_similarity_comparison.py`](visualize_similarity_comparison.py) visualizes comparison of feature similarity between target pairs with common functions and those without common functions.

+ [`parse_gtex_tpm.R`](parse_gtex_tpm.R) processes raw GTEx tissue-specific expression data for gene ID mapping and differential expression analysis. 

+ [`analyze_select_target_expression.R`](analyze_select_target_expression.R) analyzes whether selected target genes of adverse events are more likely to be differentially expressed in the matched tissues using GTEx expression data.

+ [`visualize_expression_comparison.py`](visualize_expression_comparison.py) visualizes comparison of differential expression values between selected target genes and all measured genes.

+ [`parse_ctd.R`](parse_ctd.R) processes raw Comparative Toxicogenomics Database (CTD) gene-disease connection data for identifying disease terms of interest.  

+ [`analyze_select_target_disease.R`](analyze_select_target_disease.R) analyzes whether selected targets are more likely to be disease-related genes using CTD gene-diseaese connection data.

+ [`analyze_select_target_pathway.R`](analyze_select_target_pathway.R) performs Gene Ontology/Reactome pathway enrichment analysis of target genes of interest.

+ [`functions.R`](functions.R) contains R functions required for other scripts in the repository.

