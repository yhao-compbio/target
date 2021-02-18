# This folder contains source code used by the repository.

## R/python scripts 

+ [`integrate_druggable_genome.R`](integrate_druggable_genome.R) integrates three different resources to build human druggable genome. 

+ [`compute_target_function_similarity.R`](compute_target_function_similarity.R) computes the similarity of functional annotations between given target pairs.

+ [`visualize_function_comparison.py`](visualize_function_comparison.py) visualizes comparison of feature similarity between target pairs with common functions and those without common functions.

+ [`parse_gtex_tpm.R`](parse_gtex_tpm.R) processes raw GTEx tissue-specific expression data for gene ID mapping and differential expression analysis. 

+ [`analyze_select_target_expression.R`](analyze_select_target_expression.R) analyzes whether selected target genes of adverse events are more likely to be differentially expressed in the matched tissues using GTEx expression data.

+ [`visualize_expression_comparison.py`](visualize_expression_comparison.py) visualizes comparison of differential expression values between selected target genes and all measured genes.

+ [`parse_ctd.R`](parse_ctd.R) processes raw Comparative Toxicogenomics Database (CTD) gene-disease connection data for identifying disease terms of interest.  

+ [`analyze_select_target_disease.R`](analyze_select_target_disease.R) analyzes whether selected targets are more likely to be disease-related genes using CTD gene-diseaese connection data.

+ [`analyze_select_target_enrichment.R`](analyze_select_target_enrichment.R) performs Gene Ontology/Reactome pathway enrichment analysis of target genes of interest.

+ [`generate_ae_target_map.R`](generate_ae_target_map.Ri) generates map between adverse event and predictive targets.  

+ [`visualize_target_map.py`](visualize_target_map.py) visualizes relationship map between adverse events and predictive target features. 

+ [`cluster_target_map.R`](cluster_target_map.R) implements pvclust package to identify significant clusters of adverse events from adverse event-predictive targets map.

+ [`functions.R`](functions.R) contains R functions required for other scripts in the repository.

