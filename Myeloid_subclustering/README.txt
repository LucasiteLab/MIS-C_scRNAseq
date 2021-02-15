# Myeloid subclustering

The myeloid UMAP and subsequent analysis (e.g. computing module scores and assessing individual gene expression) can be done once the subclustering Seurat object is made. 

pbmc_cluster_annotation.csv provides the cluster number and annotation of each of our clusters in the PBMC seurat object, which can be used for subclustering and further annotation.
Serum_enrichment computes an enrichment of myeloid annotated genes in the serum data from the differential expression analysis of MIS-C vs C.HD serum. 
