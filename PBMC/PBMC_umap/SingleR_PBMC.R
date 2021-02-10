library(Seurat)
library(ggplot2)
library(SingleR)

# Single R

misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_updated_cond.rds")

monaco <- MonacoImmuneData(ensembl=FALSE)

# convert to sce and run SingleR to annotate specific clusters

misc.sce <- as.SingleCellExperiment(misc.cluster, assay = 'integrated')
annot_clust <- SingleR(test = misc.sce, ref = monaco, 
                       labels = monaco$label.fine)

misc.cluster[["SingleR.labels"]] <- annot_clust$labels

saveRDS(misc.cluster, file = "RC_SingleR/misc.cluster_singleR.rds")


# UMAP
plot8 <- DimPlot(misc.cluster, reduction = "umap", 
                 cols=sample(DiscretePalette(30, palette='polychrome')),
                 pt.size =0.05,
                 group.by = "SingleR.labels", label = TRUE)
ggsave("RC_SingleR/pbmc_labels_dimplot_graphlabels.pdf", plot = plot8, height =12, width = 16)