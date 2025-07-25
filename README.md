RNAseq.R: A Comprehensive Pipeline for RNA Sequencing Differential Expression Analysis


Introduction
RNAseq.R is an R script designed for the analysis of RNA sequencing (RNA-seq) data to identify differentially expressed genes across various experimental conditions. Leveraging the Bioconductor ecosystem, particularly the edgeR package, this script performs a full workflow including data preprocessing, differential expression analysis, visualization, and downstream functional enrichment analyses such as Gene Ontology (GO) and KEGG pathway analysis. The script is tailored for mouse (Mus musculus) RNA-seq data, as demonstrated by its use of the org.Mm.eg.db annotation package and the GSE60450 dataset.
This tool is intended for bioinformaticians and computational biologists with a strong background in R programming and statistical analysis of high-throughput sequencing data. It provides a robust framework for analyzing count-based RNA-seq data, generating publication-quality visualizations, and performing advanced gene set enrichment analyses.

Prerequisites
Software Requirements

R: Version 4.0.0 or higher recommended.
Bioconductor: Version 3.12 or higher.

Required R Packages
The script depends on the following R packages, which should be installed prior to execution:

Core analysis packages:
BiocManager
edgeR
limma
RnaSeqGeneEdgeRQL
org.Mm.eg.db (for gene annotation)


Visualization packages:
ggplot2
pheatmap
heatmaply
RColorBrewer


Functional enrichment packages:
clusterProfiler
enrichplot
GO.db
KEGGREST
UpSetR


Network analysis:
WGCNA
impute
preprocessCore
dynamicTreeCut
fastcluster


Miscellaneous:
dplyr
tidyr
RCurl



To install these dependencies, run:
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "RnaSeqGeneEdgeRQL", "org.Mm.eg.db", "clusterProfiler", "enrichplot", "GO.db", "KEGGREST", "WGCNA", "impute", "preprocessCore", "dynamicTreeCut", "fastcluster", "ggplot2", "pheatmap", "heatmaply", "RColorBrewer", "UpSetR", "dplyr", "tidyr", "RCurl"))


Input Data
Data Format
The script expects RNA-seq count data in a tabular format:

Rows: Genes, identified by Entrez Gene IDs.
Columns: Samples, with sample identifiers as column names.
Content: Raw read counts per gene per sample.

Example Dataset
The script downloads and uses the GSE60450 dataset (GSE60450_Lactation-GenewiseCounts.txt.gz) from GEO if not already present locally. This dataset contains gene-wise count data for mouse mammary gland samples across different physiological states (virgin, pregnant, lactating) and cell types (basal and luminal).
Targets File
A targets file (targets.txt) is required to define the experimental design:

Columns: Typically include CellType (e.g., B for basal, L for luminal) and Status (e.g., virgin, pregnant, lactating).
Format: Tab-delimited text file.
Source: Provided via the RnaSeqGeneEdgeRQL package (system.file("extdata", "targets.txt", package="RnaSeqGeneEdgeRQL")).


Script Overview
The RNAseq.R script executes the following major steps:

Data Loading and Preprocessing: Loads count data, annotates genes, filters lowly expressed genes, and normalizes the data.
Differential Expression Analysis: Uses edgeR to perform statistical testing for differential expression across specified contrasts.
Visualization: Generates diagnostic plots (e.g., MDS, MD plots, volcano plots) and heatmaps.
Functional Enrichment: Conducts GO and KEGG pathway analyses, including over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
Network Analysis: Performs Weighted Gene Co-expression Network Analysis (WGCNA) to identify co-expression modules.
Output Generation: Saves results and plots to structured directories.


Detailed Steps
1. Data Loading and Preprocessing

Loading: Reads count data using read.delim and the targets file to define sample groups (e.g., B.lactating, L.virgin).
Annotation: Maps Entrez Gene IDs to gene symbols using org.Mm.eg.db.
Filtering: Removes genes with missing symbols and applies filterByExpr to retain only genes with sufficient expression levels.
Normalization: Calculates normalization factors using the Trimmed Mean of M-values (TMM) method via calcNormFactors.

2. Differential Expression Analysis

Design Matrix: Constructs a model matrix based on sample groups (e.g., ~0+group).
Dispersion Estimation: Estimates dispersion using estimateDisp with robust settings.
Model Fitting: Fits a quasi-likelihood negative binomial model using glmQLFit.
Contrasts: Defines pairwise contrasts (e.g., B.lactating - B.pregnant) and tests for differential expression using glmQLFTest and glmTreat (for log-fold change thresholding).
Advanced Tests: Includes ANOVA-like tests and differential-of-differences analyses.

3. Visualization

MDS Plot: Visualizes sample relationships using plotMDS.
MD Plots: Displays mean-difference plots for each sample and contrast.
Volcano Plots: Highlights significant differentially expressed genes per contrast.
Heatmaps: Generates static (pheatmap) and interactive (heatmaply) heatmaps for top differentially expressed genes.

4. Functional Enrichment

GO Analysis: Performs ORA (goana) and GSEA (gseGO) using clusterProfiler, focusing on Biological Process (BP) terms.
KEGG Analysis: Conducts ORA (enrichKEGG) and GSEA (gseKEGG) for pathway enrichment.
Visualization: Produces ridge plots, dot plots, and detailed GSEA plots for significant terms/pathways.

5. Network Analysis

WGCNA: Identifies co-expression modules by:
Filtering genes based on variance.
Selecting a soft-thresholding power using pickSoftThreshold.
Generating diagnostic plots for power selection.



6. Output Generation
Results are saved in a structured directory hierarchy (e.g., volcano_plots/, GO_analysis/) as CSV files and plots (PNG/PDF).

Output
The script generates the following output files and directories:
CSV Files

gene_counts.csv: Raw count data.
normalized_samples.csv: Normalization factors and library sizes.
design_matrix.csv: Experimental design matrix.
filtered_expression_data.csv: Filtered and normalized count data.
glm_coefficients.csv: Coefficients from the fitted model.
all_DE_results.csv: Comprehensive differential expression results.
[contrast]_DE_results.csv: Results for specific contrasts (e.g., B_LvsP_DE_results.csv).
Enrichment results: Various CSV files in GO_analysis/, KEGG_analysis/, etc.

Plot Directories

volcano_plots/: Volcano plots for each contrast.
md_plots_annotated/: Annotated MD plots.
treat_volcano/: Volcano plots for glmTreat results.
md_treat_plots/: MD plots for glmTreat results.
heatmap_all_contrasts/: Heatmaps across contrasts.
GO_analysis/: GO enrichment plots (ridge, dot, detailed).
KEGG_analysis/: KEGG enrichment plots and results.
WGCNA_results/: WGCNA power selection plots.


Usage
To execute the script:

Ensure all dependencies are installed.
Place input files (if not using GSE60450) in the working directory or modify file paths accordingly.
Run the script in an R environment:source("RNAseq.R")



No command-line arguments are required, as the script uses hardcoded file paths and parameters. Modify the script directly to adjust contrasts, thresholds, or file locations.

Troubleshooting

Package Installation Errors: Use BiocManager::install() to resolve missing packages. Check R and Bioconductor versions.
File Not Found: Verify the GSE60450 file is downloaded or provide a local copy. Adjust paths in the script if necessary.
Lowly Expressed Genes Warning: This is expected behavior during filtering and does not indicate an error.
Memory Issues: For large datasets, increase R’s memory allocation or filter more aggressively.


References

Chen Y, Lun ATL, Smyth GK (2016). “From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.” F1000Research, 5, 1438. http://f1000research.com/articles/5-1438.
Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140. DOI: 10.1093/bioinformatics/btp616
Ritchie, M. E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. DOI: 10.1093/nar/gkv007
Wu, T., et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2(3), 100141. DOI: 10.1016/j.xinn.2021.100141
Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559. DOI: 10.1186/1471-2105-9-559
