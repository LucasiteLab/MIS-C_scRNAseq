# T and NK subclustering

The T and NK UMAP and subsequent analysis can be done starting with the PBMC Seurat object. The pbmc_cluster_annotation.csv (main) provides the cluster number and annotation of each of our clusters in the PBMC seurat object. 

- ADT UMAP overlay: see A1 (umap_ADT_GEX_overlay_plots.R) ** Fig S4B
- Dot plots for cluster annotation: see A1 (Dotplot_cluster_markers.R) ** Fig 4B, S4A

T.NK_umap.R ** Fig 4A
T.NK_cytotoxicity ** Fig 4E-F, S4G, S4H
FindMarkers_heatmaps.R ** Fig 4D, S4D (CD8 memory subset used as representative example)
Cell_type_distribution_T.NK.R ** Fig 4C


# TCR analysis

TCR preprocessing follows the workflow shown in TCR_analysis/Preprocessing_workflow. As in GEX, TCR data is incorporated from different studies and pre-processed in parallel. The input for this data is filtered_contig_annotation.csv data found in the GSE associated with our manuscript. 

Abundance_Diversity_Analysis.R ** Fig S4I (input is the output from TCR preprocessing. diversity.csv is output)
Diversity_plot_non_naive.R ** Fig S4I (diversity.csv is input)


