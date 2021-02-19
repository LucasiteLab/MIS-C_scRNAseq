# Myeloid subclustering

The myeloid UMAP and subsequent analysis (e.g. computing module scores and assessing individual gene expression) is done based on the PBMC Seurat object (See A1)
pbmc_cluster_annotation.csv (in main page) provides the cluster number and annotation of each of our clusters in the PBMC Seurat object. 
Serum_enrichment computes an enrichment of myeloid annotated genes in the serum data from the differential expression analysis of MIS-C vs C.HD serum. 

-A representative DEG analysis and heatmap construction can be found in A3 (Findmarkers_heatmaps.R) ** Fig 3C
-Dotplots for cluster annotation can be found in A1 (Dotplot_cluster_markers.R) ** Fig 3B
-Pathway enrichment analyses were done using the online web tool EnrichR


Myeloid_umap.R ** Fig 3A
Myeloid_alarmin_HLA_CD86.R ** Fig D-E, Fig S3B-D
ModuleScore_sepsis.R ** Fig 3G
Serum_volcano_enrichment.R ** Fig 3H
Cell_type_distribution_myeloid.R ** Fig S3A