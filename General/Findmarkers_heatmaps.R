library(Seurat)
library(ggplot2)
library(tidyverse)


misc.cluster <- readRDS("shared/MISC_NKT.integrated_090620.rds") # any seurat object
Idents(misc.cluster) <- "subcluster_annotation1"


# Effector memory/memory subsets
cells.to.keep <- WhichCells(misc.cluster, idents = c("Terminal_effector_memory_CD8_Tcells", "Central_memory_CD8_Tcells", 
                                                     "Effector_memory_CD8_Tcells") ) 
misc.mem <- subset(misc.cluster, cells = cells.to.keep)
markers.mem <- FindMarkers(misc.mem, group.by="condition", ident.1="MIS-C", ident.2="C.HD", assay = "RNA")
write.csv(markers.mem, file = "markers_combined_CD8_effector_memory.csv")


markers_mem_up <- markers_mem %>% filter(avg_logFC > 0.5  & p_val_adj < 0.05)
markers_mem_dn <- markers_mem %>% filter(avg_logFC < -0.5  & p_val_adj < 0.05)
markers_mem_up <- markers_mem_up %>% arrange(desc(avg_logFC))
markers_mem_dn <- markers_mem_dn %>% arrange(avg_logFC)
top.markers_mem_up <- head(markers_mem_up, 20)
top.markers_mem_dn <- head(markers_mem_dn, 20)


# T-cell
misc.tcell <- readRDS("shared/t_cell/misc_integrated_tcell_updated_cond.rds")
misc.cd8 <- subset(misc.tcell, idents= c("Effector_memory_CD8_Tcells", 
                                         "Terminal_effector_memory_CD8_Tcells", 
                                         "Central_memory_CD8_Tcells"))


#CD8
features.to.use <- c(top.markers_mem_up$X, top.markers_mem_dn$X)
Idents(misc.cd8) <- "sample_id"
avg.topmarkers <- AverageExpression(misc.cd8, assays = "RNA", 
                                    #slot= "scale.data", 
                                    features = features.to.use)

col_order <-  c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6",
                "P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1",
                "P3.2", "P4.2", "A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                "A.HD13", "A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1", "A.COV1.2",
                "A.COV2.2", "A.COV3.2", "A.COV4.2", 
                "A.COV5.2", 
                "A.COV6.2")

avg.topmarkers$RNA <- avg.topmarkers$RNA[, col_order]
avg.topmarkers$RNA[16:ncol(avg.topmarkers$RNA)] <- NULL
avg.markers1 <- data.matrix(avg.topmarkers[[1]])

avg.markers1 <- genescale(avg.markers1, axis=1, method="Z")

# annotations
col_groups1 <- c(rep("C.HD",6), rep("MIS-C",7), rep("MIS-C-R",2))

annotations <- data.frame(group = col_groups1)
rownames(annotations) <- colnames(avg.topmarkers$RNA)
annotations_t <- t(as.matrix(annotations))
mat_colors <- list(group =  cols <- c("#6baed6", "#c94040", "#969696"))
names(mat_colors$group) <- unique(col_groups1)


pdf(file = "paper_heatmaps/complexHM_cd8_2.0.pdf", height = 8, width = 3)

Heatmap(avg.markers, col = inferno(10), row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"), 
        column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica")) %v% 
  Heatmap(avg.markers1, col = inferno(10),  row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"), 
          column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"))  %v% 
  Heatmap(annotations_t, col = mat_colors$group, column_names_gp = gpar(fontsize =5, fontfamily = "Helvetica")) 

dev.off()