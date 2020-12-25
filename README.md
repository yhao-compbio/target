# The aim of this repository is to analyze the biological connections between targets and drug toxicity. 

## Introduction 

In the [TTox repository](https://github.com/yhao-compbio/TTox), we developed a feature selection pipeline that can identify target proteins predictive of drug toxicity. An advantage of our pipeline is that the identified targets can help explain explain cellular mechanisms underlying structure-toxicity associations. In this analysis, we uncovered the connections between identified targets and drug toxicity, thus demonstrated the significance of our feature selection results from a biological perspective. 

The codes to perform our analyses can be found in [`src/`](src/), along with detailed documentation.

## Results 

### Our pipeline can identify similar structure properties for target proteins of similar function

First, we evaluated the predictive molecular descriptors that our pipeline identified for target binding. Given there is no benchmarking standards connection descriptors to targets, we evaluated the results by checking whether similar descriptors are identified for targets of similar function. We grouped 367 target proteins (with model AUROC > 0.85) into six function classes: enzyme (33.8%), G-protein coupled receptor (24.3%), catalytic receptor (18.5%), nuclear hormone receptor (4.9%), ion channel (4.1%), transporter (3.5%). The pie chart below (on the left) shows the distribution of function classes. We computed pairwise Jaccard similarity of identified descriptors among targets, compared the similarity within and beyond each class. The boxplot below (on the right) shows the comparison. For G-protein coupled receptors and catalytic receptors, the intraclass similarity is significantly higher than the interclass similarity. Note that these two categories account for 44 percent of the targets we analyzed. 

![class_comparison](https://github.com/yhao-compbio/TTox/blob/master/plot/compound_target_0.25_binary_feature_select_implementation/descriptor_all_compare/class_comparison_descriptor_all.png)

We also mapped the targets to their Gene Ontology and Reactome pathway annotations. Detailed documentation about the annotations can be found at the [ontology repository](https://github.com/yhao-compbio/ontology). Again, we compared the similarity between target pairs with common function to those without common function. The boxplot below shows the comparison. We found that target pairs with common function are more likely to be identified with similar descriptors.  

![function_similarity_group](plot/target_similarity/descriptor_all_select_features_mc_0.85_function_similarity_group_boxplot.png)

### Identified target genes are differentially expressed in the tissue of toxicity 

Next, we studied the differential expression of target genes in the tissue of toxicity. For this analysis, we used GTEx tissue-specific TPM data. We normalized the expression values of each gene by the median expression across all tissues, then took the absolute value of logarithms as a measure of differential expression. We mapped 12 adverse events of interest (with model AUROC > 0.65) to a GTEx tissue, compared the differential expression of identified target genes to the null distribution consisting of all genes. The boxplot below shows the comparison. The identified target genes of eight adverse events are differentially expressed in their tissues of toxicity. They are: cardiac arrest (heart), myoclonus (muscle), pseudomembranous colitis (colon), nephropathy (kidney), renal tubular disorder (kidney), rash erythematous (skin), rash macular (skin), and respiratory acidosis (lung). 

![select_targets_de_compare](plot/target_expression/descriptor_all_all_adverse_event_select_targets_de_compare.png)

## Identified target genes are enriched for disease markers/therapeutics

Last, we studied the connections between target genes and associated diseases. For this analysis, we used disease-gene connections from CTD. We only kept the connections with direct supporting evidence (e.g. the gene is a disease marker/therapeutic). We mapped 12 adverse events of interest (with model AUROC > 0.65) to a list of associated diseases through key words matching. We then tested whether the identified target genes are enriched for disease-related genes. The table below shows the results of Fisher's exact test. In eight of the 12 adverse events, the identified target genes are twice or more likely to be disease-related genes (OR > 2). In three of the 12 adverse events, the identified target genes are significantly enriched for disease-related genes (FDR < 0.05). They are: arrhythmia-heart disease, rash erythematous-skin disease, and rash macular-skin disease. 

![select_targets_disease](data/target_disease/descriptor_all_all_adverse_event_select_features_symbol_connection_enrich.png)

Detailed documentation about our results and figures can be found at [`data/`](data/) and [`plot/`](plot/), respectively.

## Conclusion
Through these analyses, we demonstrated that not only can our feature selection pipeline identify features predictive of the outcome, it can also identify features with biological significance. 

## References

+ Kumar RD, Chang LW, Ellis MJ, Bose R. Prioritizing potentially druggable mutations with dGene: an annotation tool for cancer genome sequencing data. PloS one. 2013 Jun 27;8(6):e67980.

+ Armstrong JF, Faccenda E, Harding SD, Pawson AJ, Southan C, Sharman JL, Campo B, Cavanagh DR, Alexander SP, Davenport AP, Spedding M. The IUPHAR/BPS Guide to PHARMACOLOGY in 2020: extending immunopharmacology content and introducing the IUPHAR/MMV Guide to MALARIA PHARMACOLOGY. Nucleic acids research. 2020 Jan 8;48(D1):D1006-21.

+ Wishart DS, Feunang YD, Guo AC, Lo EJ, Marcu A, Grant JR, Sajed T, Johnson D, Li C, Sayeeda Z, Assempour N. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic acids research. 2018 Jan 4;46(D1):D1074-82
