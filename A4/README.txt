# B subclustering


The B cell UMAP and subsequent analysis (e.g.assessing gene expression or BCR analysis with annotated clusters) is done starting from PBMC Seurat object. The pbmc_cluster_annotation.csv (main) provides the cluster number and annotation of each of our clusters in the PBMC seurat object. 

- Dot plots for cluster annotation as in A1 (Dotplot_cluster_markers.R)

Bcell_umap.R ** Fig 5A
T-B_cell_correlation.R (for both scRNA-seq and flow- data included) ** Fig 5G
Cell_type_distribution_bcell.R ** Fig 5C


# BCR Analysis

Scripts for performing BCR sequence analysis
** Fig 5D-F, S5F-G, 6D-E, S6F

Requirements:
Immcantation 4.0.0 Docker image
Change-O v1.0.1

R packages:
alakazam v1.1.0.999
shazam v1.0.2.999
dplyr v1.0.2
tidyr v1.1.2
stringr v1.4.0
stringdist v0.9.6.3
ggplot2 v3.3.2 
ggpubr v0.4.0

To reproduce analyses:
-Run the Immcantation docker image
-Run processData.sh to align BCR sequences to IMGT reference and convert to AIRR format.
-Run combineData.R to compile data into a single table with appropriate metadata
-Run cloneGermline.R to group BCRs into clonal clusters and reconstruct clonal germlines
-Run analysis.R to perform analyses and generate figures*.

* Relevant Figures for BCR analysis are relevant for both A4 and A5
