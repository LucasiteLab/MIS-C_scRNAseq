Anjali Ramaswamy
Lucas Lab
2/10/2021


Scripts for performing PBMC clustering, sub-clustering of myeloid, T.NK cells, and B cells. Also included are scripts for:
GEX and ADT:
- SingleR and ADT (antibody dependent tag) UMAP overlay.
- Performing FindMarkers DEG analysis and heatmap generation.
- ModuleScore signature scoring and AUCell signature enrichment for viral/bacterial signatures, alarmin signatures, and similar. 
- Plotting specific genes and signatures, e.g. scaled average expression of cytotoxicity markers.
- EBV and CMV alignment.
- Connectome analysis. 

Antigen receptor analysis:
- Performing BCR and TCR analyses 

Relevant data and instructions are found in each subdirectory, and pre-processed data for pediatric samples and A.HD1-3 can be found at GSE166489 (e.g. filtered barcodes, genes, and matrix files, along with TCR/BCR filtered contig annotation files for each sample. 

A1 -- Figure 2, S2
A2 -- Figure 3, S3
A3 -- Figure 4, S4
A4 -- Figure 5, S5
A5 -- Figure 6, S6

R packages:

Seurat v3.2.1 
AUCell v1.12.0
ggplot2 v3.3.2
SingleR v1.4.0
ComplexHeatmap v2.5.5 
stats v4.0.2
Connectome v0.2.2 
alakazam v1.1.0.999
shazam v1.0.2.999
dplyr v1.0.2
tidyr v1.1.2
stringr v1.4.0
stringdist v0.9.6.3
ggpubr v0.4.0


