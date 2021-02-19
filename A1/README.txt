A1

Input needed are hashing_files and filtered_feature_bc_matrix.tar.gz (These processed files are stored in GSE166489). 
A.HD4-A.HD7 along with all adult COVID-19 samples can be found corresponding to Unterman et al. medRxiv 2020. A.HD8-A.HD13 can be found corresponding to Pappalardo et al. Science Immunology 2020. 

# PBMC clustering

PBMC_umap.R; SingleR_PBMC.R; PBMC_umap_ADT_layer.R
To generate PBMC UMAP, use filtered_feature_bc_matrix.tar.gz files as input for each sample, as well as hashing files (relevant for hashed COVID-19 and Adult healthy donor samples). Datasets are filtered, integrated, and clustered. Dead cell clusters are identified and removed, and the data are re-integrated. Also included is code for the SingleR automatic annotation, and adding a CITE-seq layer. These tools were used in conjunction with cluster differential expressed genes to inform cluster annotations. 
** Fig 2A


# Subsequent analysis
Dotplot_cluster_markers.R ** Fig 2B, S2A
Cell_type_distribution_grid.R ** Fig 2C
ModuleScore_viral_bacterial.R; AUCell_type_I_ifn.R ** Fig 2D, S2G
umap_ADT_GEX_overlay_plots.R ** Fig S2B, S2C
Connectome_pbmc.R ** Fig S2F
chk_viral_counts.py; ebv_cmv.sh ** Fig S2H
Correlation_plot.R ** Fig S2D







