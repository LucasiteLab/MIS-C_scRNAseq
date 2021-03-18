A1

Input needed are hashing_files and filtered_feature_bc_matrix.tar.gz (These processed files are stored in GSE166489). 
A.HD4-A.HD7 along with all adult COVID-19 samples can be found corresponding to Unterman et al. medRxiv 2020. A.HD8-A.HD13 can be found corresponding to Pappalardo et al. Science Immunology 2020. 

# PBMC clustering

To generate PBMC UMAP (PBMC_umap.R), use filtered_feature_bc_matrix.tar.gz files as input for each sample, as well as hashing files (relevant for hashed COVID-19 and Adult healthy donor samples). Datasets are filtered, integrated, and clustered. Dead cell clusters are identified and removed, and the data are re-integrated. Also included is code for the SingleR automatic annotation, and adding a CITE-seq layer. These tools were used in conjunction with cluster differential expressed genes to inform cluster annotations. 

PBMC_umap.R; SingleR_PBMC.R;  Fig 2A starting from CellRanger output
Ramaswamy2021_MIS-C_10x_Bcell_3182021.ipynb ** Fig 2A starting from processed objects (FASTgenomics)
PBMC_umap_ADT_layer.R ** Generating data layer for S2B-C


# Subsequent analysis
Ramaswamy2021_MIS-C_10x_Bcell_3182021.ipynb ** Fig 2B-D, S2A
AUCell_type_I_ifn.R ** Fig S2G
umap_ADT_GEX_overlay_plots.R ** Fig S2B, S2C
Connectome_pbmc.R ** Fig S2F
chk_viral_counts.py; ebv_cmv.sh ** Fig S2H
Correlation_plot.R ** Fig S2D







