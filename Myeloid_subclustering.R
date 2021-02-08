library(dplyr)
library(Seurat)
library(cowplot)
library(future)
library(ggplot2)
library(readr)
library(patchwork)
library(SingleR)
library(reshape2)
library(tidyverse)

misc.integrated.new  <- readRDS(file="~/misc_integrated_object_2.0.rds")

#Subset for platelets, pDC and myeloid cells:
myeloid <- subset(misc.integrated.new, idents = c(2,8,12,15,21,23,24,25,27))
Idents(myeloid) <- "orig.ident"
orig.composition <- table(Idents (object = myeloid), myeloid@meta.data$annotation)
write.table(orig.composition, file="orig.composition.txt", quote=F, sep="\t", na="", row.names=T, col.names=T)

#Re-integrate the object, with sample 10 and 14 as reference:
DefaultAssay(myeloid) <- "RNA"
myeloid[["integrated"]] <- NULL
myeloid.list <- SplitObject(myeloid, split.by = "orig.ident")
myeloid.list <- myeloid.list[c("NS1A","NS1B","TS2A","TS2B","TS3A","TS3B","TP9B","NS0A", "NS0B","Y111-1","Y113-1", "HD_35F","HD_32M","HD_36M","NC-13F","Y117-1", "Y117-R","Y124-1","Y124-R","Y125-1","Y127-1","Y129-1","Y28-2", "Y28-4","Y29-2","Y54-4", "Y70-4","HA5876","HA5877","HA5894", "HA5952","HA5957","HA5953","C39","C32","C27","C33","TP8B")]  
for (i in 1:length(myeloid.list)) {
  myeloid.list[[i]] <- NormalizeData(myeloid.list[[i]], verbose = FALSE)
  myeloid.list[[i]] <- FindVariableFeatures(myeloid.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}
myeloid.anchors <- FindIntegrationAnchors(object.list = myeloid.list, dims = 1:30, reference =c(10,14))
myeloid.integrated <- IntegrateData(anchorset = myeloid.anchors, dims = 1:30)
DefaultAssay(myeloid.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering:
myeloid.integrated <- ScaleData(object = myeloid.integrated)
myeloid.integrated <- RunPCA(object = myeloid.integrated)
ElbowPlot(object = myeloid.integrated, ndims = 25)
myeloid.integrated <- FindNeighbors(object = myeloid.integrated, dims = 1:20)
myeloid.integrated <- FindClusters(object = myeloid.integrated, resolution = 0.5)
myeloid.integrated <- RunUMAP(object = myeloid.integrated, dims = 1:20)

#Definition of cell type per cluster:
DefaultAssay(myeloid.integrated) <- "RNA"
myeloid.integrated <- ScaleData(object = myeloid.integrated)
myeloid.integrated.markers <- FindAllMarkers(object = myeloid.integrated, min.pct = 0.25, logfc.threshold = 0.25)
myeloid.integrated.markers_roc <- FindAllMarkers(object = myeloid.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
VlnPlot(object = myeloid.integrated, features = c('FCGR3A', 'CD14', 'S100A9', 'TUBB1', 'PPBP', 'PF4'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('ZBTB46', 'CD1C', 'ITGAX', 'CLEC4C','LILRA4','FCER1A'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('CD4','CD8A','CD3E','IL7R', 'CCR7', 'S100A4'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('FCGR3A','NCAM1', 'GNLY', 'NKG7', 'CD3E', 'CD8A'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('CD4','CD8A','CD3E','IL7R', 'CCR7', 'S100A4'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('FCGR3A','NCAM1','GNLY', 'KLRB1', 'PRF1', 'GZMA'), pt.size = 0.2)
VlnPlot(object = myeloid.integrated, features = c('MS4A1','CD79A', 'IGHG2', 'CD38', 'CD27', 'IGHD'), pt.size = 0.2)

#Filter out contaminants and final visualization:
DefaultAssay(myeloid.integrated) <- "RNA"
myeloid.integrated.filtered <- subset(x = myeloid.integrated, idents = c(0,1,2,3,4,7,9,11,12,15))
DimPlot(object = myeloid.integrated.filtered, reduction = 'umap', label = TRUE)
new.cluster.ids <- c("Neutrophils I", "Neutrophils II","Classical Monocytes I","Classical Monocytes II", "Non Classical Monocytes", "Intermediate Monocytes", "Conventional DC", "Plasmacytoid DC", "ISG+ Myeloid Cells", "Bad Quality Cells")
names(x = new.cluster.ids) <- levels(x = myeloid.integrated.filtered)
myeloid.integrated.filtered <- RenameIdents(object = myeloid.integrated.filtered, new.cluster.ids)
DimPlot(object = myeloid.integrated.filtered, reduction = 'umap', cols = c("chartreuse3", "springgreen4",  "plum1", "maroon1", "purple",   "deeppink3","dodgerblue", "deepskyblue1","orchid1",  "dimgrey"))
table(Idents (object = myeloid.integrated.filtered), myeloid.integrated.filtered@meta.data$condition)
table(Idents (object = myeloid.integrated.filtered), myeloid.integrated.filtered@meta.data$orig.ident)
table(Idents (object = myeloid.integrated.filtered), myeloid.integrated.filtered@meta.data$storage)

#Validation with SingleR package:
test <- GetAssayData(myeloid.integrated.filtered)
monaco <- MonacoImmuneData(ensembl=FALSE)
#For single cells
annot <- SingleR(test = test, ref = monaco, labels = monaco$label.fine)
myeloid.integrated.filtered[["SingleR.labels"]] <- annot$labels
DimPlot(myeloid.integrated.filtered, reduction = "umap",pt.size =0.05, group.by = "SingleR.labels", cols = c("maroon1", "purple",  "springgreen4",  "dodgerblue","deepskyblue1"))
#For defined clusters
myeloid.integrated.filtered.sce <- as.SingleCellExperiment(myeloid.integrated.filtered, assay = 'integrated')
annot_clust <- SingleR(test = myeloid.integrated.filtered.sce, ref = monaco, labels = monaco$label.fine, method = "cluster",clusters=myeloid.integrated.filtered.sce$seurat_clusters)
myeloid.integrated.filtered[["SingleR.cluster.labels"]] <- annot_clust$labels[match(myeloid.integrated.filtered[[]][["seurat_clusters"]], rownames(annot_clust))]
DimPlot(myeloid.integrated.filtered, reduction = "umap",pt.size =0.05,group.by = "SingleR.cluster.labels", cols = c("maroon1", "purple",  "springgreen4",  "dodgerblue","deepskyblue1"))