library(Seurat)
library(ggplot2)
library(viridis)

## PBMC

misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_updated_cond.rds")

exhaustive.pbmc <- c("rna_CD3D", "rna_TCF7", "rna_CD4", "rna_CCR7", "rna_FOXP3",
                     "rna_CXCR3", "rna_CD8A", "rna_MKI67", "rna_PPBP", "rna_TRAV1-2", 
                     "rna_TRDV2", "rna_GZMA", "rna_NCAM1", "rna_FCGR3A", "rna_MS4A1", 
                     "rna_IGHD", "rna_CD27", "rna_SDC1", "rna_S100A8", "rna_CD14", 
                     "rna_HLA-DRA", "rna_CD1C", "rna_CLEC4C", "rna_CD34")

features.pbmc.small <- c("rna_CD3D", "rna_TCF7", "rna_CD8A", "rna_GZMA", "rna_MKI67", "rna_PPBP", 
                         "rna_TRDV2", "rna_NCAM1","rna_IGHD", "rna_MS4A1", "rna_SDC1", "rna_CD14",
                         "rna_S100A8", "rna_FCGR3A","rna_CD1C", "rna_CD34")

level.order <- c("CD4 naive I", "CD4 naive II", "CD4 and CD8 mixed naive", "CD8 naive", "CD4 memory", 
                 "Regulatory T cells", "CD4 and CD8 activated memory", "CD8 memory", "CD8 effector", 
                 "Proliferating T and NK cell", "NK-T doublets", "T-NK-B cell doublets", 
                 "Platelet-T cell doublets", "Platelets", "Platelet-bound monocytes",
                 "MAIT-NKT cells", "gdT cells", "CD56dim CD16bright NK", "CD56bright CD16dim NK",
                 "NK-monocyte doublets", "Naive B", "Memory B", "Plasma cells", "Neutrophils", 
                 "Classical monocytes", "Non classical monocytes", "conventional DC", 
                 "plasmacytoid DC", "Erythroid ", "Progenitor")

misc.cluster$annotation <- factor(misc.cluster$annotation, levels = level.order )

Idents(misc.cluster) <- "annotation"

inferno_mod <- inferno(20)[3:18]

dot <- DotPlot(misc.cluster, features = rev(exhaustive.pbmc),  
               assay = "RNA", dot.scale = 3.75)+ 
  scale_color_gradient2(low = inferno_mod[3], mid = "#D64B40FF", high = "#F1ED6FFF", midpoint = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
        axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
  coord_flip()

ggsave(plot = dot, file = "annotation_dotplot/pbmc_exhaustive_dp_inferno.pdf", height = 6, width = 7.5)


## Myeloid

misc.cluster <- readRDS("shared/Myeloid/misc_integrated_myeloid_updated_cond.rds")

features.myeloid <- c("rna_S100A8", "rna_CSF3R", "rna_FPR1", "rna_S100A9", 
                      "rna_CD14", "rna_VCAN", "rna_CD163", "rna_FCGR3A",
                      "rna_HLA-DRA", "rna_CD1C", "rna_LILRA4")

Idents(misc.cluster) <- "FinalAnnotation"

inferno_mod <- inferno(20)[3:18]

dot <- DotPlot(misc.cluster, features = rev(features.myeloid),  
               assay = "RNA", dot.scale = 4.5)+ 
  scale_color_gradient2(low = inferno_mod[3], mid = "#D64B40FF", high = "#F1ED6FFF", midpoint = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
        axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
  coord_flip()


ggsave(plot = dot, file = "annotation_dotplot/myeloid_dp_inferno.pdf", height = 5, width = 5)

## T-cell 

misc.cluster <- readRDS("shared/t_cell/misc_integrated_tcell_updated_cond.rds")

level.order <- c("Naive_CD4_Tcells_1", "Naive_CD4_Tcells_2", "Naive_CD4_Tcells_3", "Naive_CD8_Tcells",
                 "Central_memory_CD8_Tcells", "Memory_CD4_Tcells", "CCR6_Memory_CD4_Tcells",
                 "CXCR3_Memory_CD4_Tcells", "regulatory_Tcells", "Effector_memory_CD8_Tcells", 
                 "Terminal_effector_memory_CD8_Tcells", "Ki67_NK/Tcells", "MAIT_cells", 
                 "Vd2_gdTcells", "CD56dim_S100A4_NKcells", "CD56dim_CD38NKcells",
                 "CD56bright_NKcells", "T/NK/monocyte_doublets")
misc.cluster$subcluster_annotation1 <- factor(misc.cluster$subcluster_annotation1, levels = level.order )

Idents(misc.cluster) <- "subcluster_annotation1"
misc.cluster <- subset(misc.cluster, idents = c("T/NK/monocyte_doublets"), invert = TRUE)


features.tcell.small <- c("rna_CD3D", "rna_IL7R", "rna_CD8A", "rna_CD4", "rna_CCR6", 
                          "rna_CXCR3", "rna_FOXP3", "rna_GZMA", 
                          "rna_MKI67",  "rna_TRAV1-2", "rna_TRDV2",  
                          "rna_S100A4", "rna_CD38","rna_NKG7", "rna_PRF1")

features.tcell <- c("rna_CD3D", "rna_IL7R", "rna_CD8A", "rna_CD4", "rna_CCR7", "rna_SELL", "rna_TCF7",
                    "rna_LEF1",  "rna_CCR6", "rna_CXCR3", "rna_FOXP3", "rna_GZMA", "rna_GZMB",
                    "rna_GZMH","rna_PRF1", "rna_TIGIT", "rna_CTLA4", "rna_LAG3",  "rna_HAVCR2", "rna_MKI67", 
                    "rna_TRAV1-2", "rna_TRDV2", "rna_TRGV9","rna_IL2RA", 
                    "rna_FCGR3A", "rna_LYZ", "rna_S100A4","rna_CD38", "rna_NKG7")

dot <- DotPlot(misc.cluster, features = rev(features.tcell.small),  
               assay = "RNA", dot.scale = 4.5)+ 
  scale_color_gradient2(low = inferno_mod[3], mid = "#D64B40FF", high = "#F1ED6FFF", midpoint = 0.25) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
        axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
  coord_flip()

ggsave(plot = dot, file = "annotation_dotplot/tcell_dp_small_inferno.pdf", height = 6, width = 6.5)


## B cell

misc.cluster <- readRDS("shared/b_cell/misc_integrated_object_2.0_bcell_final.rds")


level.order <- c("Naive B I", "Naive B II", "Intermediate memory", "CD27+IgM+IgG+ memory",
                 "Memory B-cell", "Non-dividing plasma cell", "Dividing plasma cell")

misc.cluster$annotation_bcell <- factor(misc.cluster$annotation_bcell, levels = level.order )

features.bcell <- c("rna_MS4A1", "rna_IGHD", "rna_IGHM", "rna_NR4A1", 
                    "rna_IGHG3", "rna_IGHG2", "rna_IGHG1",  "rna_CD27",
                    "rna_CD38", "rna_JCHAIN", "rna_MZB1", "rna_XBP1", "rna_MKI67")

Idents(misc.cluster) <- "annotation_bcell"

dot <- DotPlot(misc.cluster, features = rev(features.bcell),  
               assay = "RNA", dot.scale = 4.5)+ 
  scale_color_gradient2(low = inferno_mod[3], mid = "#D64B40FF", high = "#F1ED6FFF", midpoint = 0.25) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
        axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
  coord_flip()

ggsave(plot = dot, file = "annotation_dotplot/bcell_dp_inferno.pdf", height = 4.5, width = 4.75)

