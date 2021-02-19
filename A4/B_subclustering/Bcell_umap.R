library(Seurat)
library(tidyverse)
library(ggplot2)
library(future)

misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0.rds") # PBMC umap
b_clusters <- WhichCells(misc.cluster, idents = c(5,10,18,22)) #Naive, mem, doublets, plasma cells
misc.bcell <- subset(misc.cluster, cells = b_clusters)

misc.bcell.list <-SplitObject(misc.bcell, split.by = "orig.ident")

sub_list <- c()

for(i in 1:length(misc.bcell.list)){
  object <- misc.bcell.list[[i]] 
  object@assays$integrated <- NULL
  DefaultAssay(object) <- "RNA"
  sub_list[[i]] <- object
}

names(sub_list) <- names(misc.bcell.list)

sub_list[[7]] <- NULL



for (i in 1:length(sub_list)) {
  sub_list[[i]] <- FindVariableFeatures(sub_list[[i]], 
                                               selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
}



misc.anchors <- FindIntegrationAnchors(object.list = sub_list, dims = 1:30, reference = c(9,13,22)) # 111_cite, 36M, 28-2 (MIS-C, A.HD, C.HD)
misc.int <- IntegrateData(anchorset = misc.anchors, dims = 1:30)

saveRDS(misc.int, file = "integrated_rds/misc_int_bcell_30dim_ref_111_28-2_36M.rds")

misc.int <- ScaleData(misc.int, verbose = FALSE)
misc.int <- RunPCA(misc.int, npcs = 30, verbose = FALSE) 
plot1 <- ElbowPlot(misc.int)
misc.cluster <- misc.int #save a copy with PCA done for UMAP and clustering later
ggsave("elbow_plots/pc_elbow_plot_int30_ref_111_28-2_36M.pdf", plot = plot1)

misc.int <- RunUMAP(misc.int, reduction = "pca", dims = 1:15) #only run the first time


plot2 <- DimPlot(misc.int, reduction = "umap", group.by = "orig.ident") #only run the first time
ggsave("umap_batch_int30_dim15_ref_111_28-2_36M.pdf", plot = plot2, height = 12, width = 14) #only run the first time


## Find neighbors and clusters (louvain) based on PCA, and then run umap 

misc.cluster <- FindNeighbors(misc.cluster, dims = 1:15) #not uMAP-ed data
misc.cluster <- FindClusters(misc.cluster, resolution = 0.3, random.seed= 10) #experiment with diff resolution
#head(Idents(bcell.cluster), 5)

colors <- c("#5A5156", "#325A9B", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
          "#B10DA1", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#E4E1E3")


misc.cluster <- RunUMAP(misc.cluster, dims = 1:15) #default is reduction = pca
plot7 <- DimPlot(misc.cluster, reduction = "umap", 
                 pt.size =0.5, label = TRUE)
ggsave("umap_louvain_res0.3int30dim15_ref_111_28-2_36M_colors.pdf", plot = plot7, height = 12, width = 14)

# Extra
toshow <- WhichCells(misc.cluster, idents = c(0))
plot8 <- DimPlot(misc.cluster, cells.highlight = toshow)
ggsave("umap_show_cluster0.png", plot = plot8)
#

saveRDS(misc.cluster, file = "misc_cov_cluster_res0.27dim15var2000_ref_111_28-2_36M_TS3B.rds")

##

bcell_markers <- FindAllMarkers(misc.cluster, only.pos = TRUE, # for annotation
                                min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

bcell_markers

write.csv(bcell_markers, file = "bcell_subcluster_markers_res0.3dim15_3ref.csv")


FindMarkers(misc.cluster, ident.1 = c(5), ident.2 = c(8), min.pct = 0.1, min.logfc = 0.1) 


## Feature plotting for annotation

plot9 <- FeaturePlot(object = misc.cluster, 
                   features = "rna_PIK3CG",
                     min.cutoff = 0,
                     #                     max.cutoff = 6,
                     pt.size = 0.8,
                     order = TRUE) + theme_classic(base_size = 25)+
  theme(axis.text = element_text(size = 35),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 45, face= "italic", hjust = 0.5)) +
  ggtitle("PIK3CG")

ggsave("pik3cg/bcell_pik3cg_umap.pdf", plot = plot9, height = 10, width = 10)



## Sample composition

cluster_prop <- table(Idents(misc.cluster), misc.cluster$orig.ident) # this goes into cell type proportion grid

write.csv(cluster_prop, file = "sample_bcell_dim15var200_3ref.csv")

## ADT

ifn_vln <- VlnPlot(misc.six, features = c("0062-anti-human-CD10"), # ADT overlay seurat object made separately
                   group.by = "seurat_clusters", 
                   pt.size = 0.5, combine = FALSE)
ifn_vln <- CombinePlots(ifn_vln)
ifn_vln <- ifn_vln + theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +ylim(c(0,400))
ggsave(ifn_vln, file = "new_bcell_umap/b_clusters_adt_CD10.png", height = 12, width = 24)


# Single R

library(SingleR)
library(ensembldb)
test <- GetAssayData(misc.cluster)
monaco <- MonacoImmuneData(ensembl=FALSE)
annot <- SingleR(test = test, ref = monaco, labels = monaco$label.fine)
misc.cluster[["SingleR.labels"]] <- annot$labels

plot8 <- DimPlot(misc.cluster, reduction = "umap", 
                 cols=DiscretePalette(32, palette='polychrome'),
                 pt.size =0.5,
                 group.by = "SingleR.labels")
ggsave("umap_bcell_louvain_labels.pdf", plot = plot8, height =12, width = 18)



saveRDS("misc_cov_cluster_res0.3dim15var2000_ref_111_28-2_36M.rds")
