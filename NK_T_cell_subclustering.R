# Download libraries
library(Seurat)
library(ggplot2)
library(data.table)
library(ggsci)

# read PBMC_files
MISC_integrated_object <- readRDS("~/misc_integrated_object_2.0.rds")

# Createseuratobject 
MISC_integrated_NKT <- subset(MISC_integrated_object, idents = c("17", "3", "20", "7", "11", "9", "13", "19", "0", "14", "16", "6", "1", "4", "26"))

Idents(MISC_integrated_NKT) <- "sample_id"


MISC_integrated_A1 <- subset(MISC_integrated_NKT, idents =  "A.COV1.1")
MISC_integrated_A12 <- subset(MISC_integrated_NKT, idents =  "A.COV1.2")
MISC_integrated_A2 <- subset(MISC_integrated_NKT, idents =  "A.COV2.1")
MISC_integrated_A22 <- subset(MISC_integrated_NKT, idents =  "A.COV2.2")
MISC_integrated_A3 <- subset(MISC_integrated_NKT, idents =  "A.COV3.1")
MISC_integrated_A32 <- subset(MISC_integrated_NKT, idents =  "A.COV3.2")
MISC_integrated_A4 <- subset(MISC_integrated_NKT, idents =  "A.COV4.1")
MISC_integrated_A42 <- subset(MISC_integrated_NKT, idents =  "A.COV4.2")
MISC_integrated_A52 <- subset(MISC_integrated_NKT, idents =  "A.COV5.2")
MISC_integrated_A62 <- subset(MISC_integrated_NKT, idents =  "A.COV6.2")
MISC_integrated_HD1 <- subset(MISC_integrated_NKT, idents =  "A.HD1")
MISC_integrated_HD2 <- subset(MISC_integrated_NKT, idents =  "A.HD2")
MISC_integrated_HD3 <- subset(MISC_integrated_NKT, idents =  "A.HD3")
MISC_integrated_HD4 <- subset(MISC_integrated_NKT, idents =  "A.HD4")
MISC_integrated_HD5 <- subset(MISC_integrated_NKT, idents =  "A.HD5")
MISC_integrated_HD6 <- subset(MISC_integrated_NKT, idents =  "A.HD6")
MISC_integrated_HD7 <- subset(MISC_integrated_NKT, idents =  "A.HD7")
MISC_integrated_HD8 <- subset(MISC_integrated_NKT, idents =  "A.HD8")
MISC_integrated_HD9 <- subset(MISC_integrated_NKT, idents =  "A.HD9")
MISC_integrated_HD10 <- subset(MISC_integrated_NKT, idents =  "A.HD10")
MISC_integrated_HD11 <- subset(MISC_integrated_NKT, idents =  "A.HD11")
MISC_integrated_HD12 <- subset(MISC_integrated_NKT, idents =  "A.HD12")
MISC_integrated_HD13 <- subset(MISC_integrated_NKT, idents =  "A.HD13")
MISC_integrated_CHD1 <- subset(MISC_integrated_NKT, idents =  "C.HD1")
MISC_integrated_CHD2 <- subset(MISC_integrated_NKT, idents =  "C.HD2")
MISC_integrated_CHD3 <- subset(MISC_integrated_NKT, idents =  "C.HD3")
MISC_integrated_CHD4 <- subset(MISC_integrated_NKT, idents =  "C.HD4")
MISC_integrated_CHD5 <- subset(MISC_integrated_NKT, idents =  "C.HD5")
MISC_integrated_CHD6 <- subset(MISC_integrated_NKT, idents =  "C.HD6")
MISC_integrated_P1 <- subset(MISC_integrated_NKT, idents =  "P1.1")
MISC_integrated_P2 <- subset(MISC_integrated_NKT, idents =  "P2.1")
MISC_integrated_P31 <- subset(MISC_integrated_NKT, idents =  "P3.1")
MISC_integrated_P32 <- subset(MISC_integrated_NKT, idents =  "P3.2")
MISC_integrated_P41 <- subset(MISC_integrated_NKT, idents =  "P4.1")
MISC_integrated_P42 <- subset(MISC_integrated_NKT, idents =  "P4.2")
MISC_integrated_P5 <- subset(MISC_integrated_NKT, idents =  "P5.1")
MISC_integrated_P6 <- subset(MISC_integrated_NKT, idents =  "P6.1")
MISC_integrated_P7 <- subset(MISC_integrated_NKT, idents =  "P7.1")

DefaultAssay(MISC_integrated_A1) <- "RNA"
DefaultAssay(MISC_integrated_A12) <- "RNA"
DefaultAssay(MISC_integrated_A2) <- "RNA"
DefaultAssay(MISC_integrated_A22) <- "RNA"
DefaultAssay(MISC_integrated_A3) <- "RNA"
DefaultAssay(MISC_integrated_A32) <- "RNA"
DefaultAssay(MISC_integrated_A4) <- "RNA"
DefaultAssay(MISC_integrated_A42) <- "RNA"
DefaultAssay(MISC_integrated_A52) <- "RNA"
DefaultAssay(MISC_integrated_A62) <- "RNA"
DefaultAssay(MISC_integrated_HD1) <- "RNA"
DefaultAssay(MISC_integrated_HD2) <- "RNA"
DefaultAssay(MISC_integrated_HD3) <- "RNA"
DefaultAssay(MISC_integrated_HD4) <- "RNA"
DefaultAssay(MISC_integrated_HD5) <- "RNA"
DefaultAssay(MISC_integrated_HD6) <- "RNA"
DefaultAssay(MISC_integrated_HD7) <- "RNA"
DefaultAssay(MISC_integrated_HD8) <- "RNA"
DefaultAssay(MISC_integrated_HD9) <- "RNA"
DefaultAssay(MISC_integrated_HD10) <- "RNA"
DefaultAssay(MISC_integrated_HD11) <- "RNA"
DefaultAssay(MISC_integrated_HD12) <- "RNA"
DefaultAssay(MISC_integrated_HD13) <- "RNA"
DefaultAssay(MISC_integrated_CHD1) <- "RNA"
DefaultAssay(MISC_integrated_CHD2) <- "RNA"
DefaultAssay(MISC_integrated_CHD3) <- "RNA"
DefaultAssay(MISC_integrated_CHD4) <- "RNA"
DefaultAssay(MISC_integrated_CHD5) <- "RNA"
DefaultAssay(MISC_integrated_CHD6) <- "RNA"
DefaultAssay(MISC_integrated_P1) <- "RNA"
DefaultAssay(MISC_integrated_P2) <- "RNA"
DefaultAssay(MISC_integrated_P31) <- "RNA"
DefaultAssay(MISC_integrated_P32) <- "RNA"
DefaultAssay(MISC_integrated_P41) <- "RNA"
DefaultAssay(MISC_integrated_P42) <- "RNA"
DefaultAssay(MISC_integrated_P5) <- "RNA"
DefaultAssay(MISC_integrated_P6) <- "RNA"
DefaultAssay(MISC_integrated_P7) <- "RNA"

MISC_integrated_A1_counts <- GetAssayData(MISC_integrated_A1, slot = "counts")
MISC_A1 <- CreateSeuratObject(counts = MISC_integrated_A1_counts, project = "MISC_A1")
MISC_integrated_A12_counts <- GetAssayData(MISC_integrated_A12, slot = "counts")
MISC_A12 <- CreateSeuratObject(counts = MISC_integrated_A12_counts, project = "MISC_A12")
MISC_integrated_A2_counts <- GetAssayData(MISC_integrated_A2, slot = "counts")
MISC_A2 <- CreateSeuratObject(counts = MISC_integrated_A2_counts, project = "MISC_A2")
MISC_integrated_A22_counts <- GetAssayData(MISC_integrated_A22, slot = "counts")
MISC_A22 <- CreateSeuratObject(counts = MISC_integrated_A22_counts, project = "MISC_A22")
MISC_integrated_A3_counts <- GetAssayData(MISC_integrated_A3, slot = "counts")
MISC_A3 <- CreateSeuratObject(counts = MISC_integrated_A3_counts, project = "MISC_A3")
MISC_integrated_A32_counts <- GetAssayData(MISC_integrated_A32, slot = "counts")
MISC_A32 <- CreateSeuratObject(counts = MISC_integrated_A32_counts, project = "MISC_A32")
MISC_integrated_A4_counts <- GetAssayData(MISC_integrated_A4, slot = "counts")
MISC_A4 <- CreateSeuratObject(counts = MISC_integrated_A4_counts, project = "MISC_A4")
MISC_integrated_A42_counts <- GetAssayData(MISC_integrated_A42, slot = "counts")
MISC_A42 <- CreateSeuratObject(counts = MISC_integrated_A42_counts, project = "MISC_A42")
MISC_integrated_A52_counts <- GetAssayData(MISC_integrated_A52, slot = "counts")
MISC_A52 <- CreateSeuratObject(counts = MISC_integrated_A52_counts, project = "MISC_A52")
MISC_integrated_A62_counts <- GetAssayData(MISC_integrated_A62, slot = "counts")
MISC_A62 <- CreateSeuratObject(counts = MISC_integrated_A62_counts, project = "MISC_A62")
MISC_integrated_HD1_counts <- GetAssayData(MISC_integrated_HD1, slot = "counts")
MISC_HD1 <- CreateSeuratObject(counts = MISC_integrated_HD1_counts, project = "MISC_HD1")
MISC_integrated_HD2_counts <- GetAssayData(MISC_integrated_HD2, slot = "counts")
MISC_HD2 <- CreateSeuratObject(counts = MISC_integrated_HD2_counts, project = "MISC_HD2")
MISC_integrated_HD3_counts <- GetAssayData(MISC_integrated_HD3, slot = "counts")
MISC_HD3 <- CreateSeuratObject(counts = MISC_integrated_HD3_counts, project = "MISC_HD3")
MISC_integrated_HD4_counts <- GetAssayData(MISC_integrated_HD4, slot = "counts")
MISC_HD4 <- CreateSeuratObject(counts = MISC_integrated_HD4_counts, project = "MISC_HD4")
MISC_integrated_HD5_counts <- GetAssayData(MISC_integrated_HD5, slot = "counts")
MISC_HD5 <- CreateSeuratObject(counts = MISC_integrated_HD5_counts, project = "MISC_HD5")
MISC_integrated_HD6_counts <- GetAssayData(MISC_integrated_HD6, slot = "counts")
MISC_HD6 <- CreateSeuratObject(counts = MISC_integrated_HD6_counts, project = "MISC_HD6")
MISC_integrated_HD7_counts <- GetAssayData(MISC_integrated_HD7, slot = "counts")
MISC_HD7 <- CreateSeuratObject(counts = MISC_integrated_HD7_counts, project = "MISC_HD7")
MISC_integrated_HD8_counts <- GetAssayData(MISC_integrated_HD8, slot = "counts")
MISC_HD8 <- CreateSeuratObject(counts = MISC_integrated_HD8_counts, project = "MISC_HD8")
MISC_integrated_HD9_counts <- GetAssayData(MISC_integrated_HD9, slot = "counts")
MISC_HD9 <- CreateSeuratObject(counts = MISC_integrated_HD9_counts, project = "MISC_HD9")
MISC_integrated_HD10_counts <- GetAssayData(MISC_integrated_HD10, slot = "counts")
MISC_HD10 <- CreateSeuratObject(counts = MISC_integrated_HD10_counts, project = "MISC_HD10")
MISC_integrated_HD11_counts <- GetAssayData(MISC_integrated_HD11, slot = "counts")
MISC_HD11 <- CreateSeuratObject(counts = MISC_integrated_HD11_counts, project = "MISC_HD11")
MISC_integrated_HD12_counts <- GetAssayData(MISC_integrated_HD12, slot = "counts")
MISC_HD12 <- CreateSeuratObject(counts = MISC_integrated_HD12_counts, project = "MISC_HD12")
MISC_integrated_HD13_counts <- GetAssayData(MISC_integrated_HD13, slot = "counts")
MISC_HD13 <- CreateSeuratObject(counts = MISC_integrated_HD13_counts, project = "MISC_HD13")
MISC_integrated_CHD1_counts <- GetAssayData(MISC_integrated_CHD1, slot = "counts")
MISC_CHD1 <- CreateSeuratObject(counts = MISC_integrated_CHD1_counts, project = "MISC_CHD1")
MISC_integrated_CHD2_counts <- GetAssayData(MISC_integrated_CHD2, slot = "counts")
MISC_CHD2 <- CreateSeuratObject(counts = MISC_integrated_CHD2_counts, project = "MISC_CHD2")
MISC_integrated_CHD3_counts <- GetAssayData(MISC_integrated_CHD3, slot = "counts")
MISC_CHD3 <- CreateSeuratObject(counts = MISC_integrated_CHD3_counts, project = "MISC_CHD3")
MISC_integrated_CHD4_counts <- GetAssayData(MISC_integrated_CHD4, slot = "counts")
MISC_CHD4 <- CreateSeuratObject(counts = MISC_integrated_CHD4_counts, project = "MISC_CHD4")
MISC_integrated_CHD5_counts <- GetAssayData(MISC_integrated_CHD5, slot = "counts")
MISC_CHD5 <- CreateSeuratObject(counts = MISC_integrated_CHD5_counts, project = "MISC_CHD5")
MISC_integrated_CHD6_counts <- GetAssayData(MISC_integrated_CHD6, slot = "counts")
MISC_CHD6 <- CreateSeuratObject(counts = MISC_integrated_CHD6_counts, project = "MISC_CHD6")
MISC_integrated_P1_counts <- GetAssayData(MISC_integrated_P1, slot = "counts")
MISC_P1 <- CreateSeuratObject(counts = MISC_integrated_P1_counts, project = "MISC_P1")
MISC_integrated_P2_counts <- GetAssayData(MISC_integrated_P2, slot = "counts")
MISC_P2 <- CreateSeuratObject(counts = MISC_integrated_P2_counts, project = "MISC_P2")
MISC_integrated_P31_counts <- GetAssayData(MISC_integrated_P31, slot = "counts")
MISC_P31 <- CreateSeuratObject(counts = MISC_integrated_P31_counts, project = "MISC_P31")
MISC_integrated_P32_counts <- GetAssayData(MISC_integrated_P32, slot = "counts")
MISC_P32 <- CreateSeuratObject(counts = MISC_integrated_P32_counts, project = "MISC_P32")
MISC_integrated_P41_counts <- GetAssayData(MISC_integrated_P41, slot = "counts")
MISC_P41 <- CreateSeuratObject(counts = MISC_integrated_P41_counts, project = "MISC_P41")
MISC_integrated_P42_counts <- GetAssayData(MISC_integrated_P42, slot = "counts")
MISC_P42 <- CreateSeuratObject(counts = MISC_integrated_P42_counts, project = "MISC_P42")
MISC_integrated_P5_counts <- GetAssayData(MISC_integrated_P5, slot = "counts")
MISC_P5 <- CreateSeuratObject(counts = MISC_integrated_P5_counts, project = "MISC_P5")
MISC_integrated_P6_counts <- GetAssayData(MISC_integrated_P6, slot = "counts")
MISC_P6 <- CreateSeuratObject(counts = MISC_integrated_P6_counts, project = "MISC_P6")
MISC_integrated_P7_counts <- GetAssayData(MISC_integrated_P7, slot = "counts")
MISC_P7 <- CreateSeuratObject(counts = MISC_integrated_P7_counts, project = "MISC_P7")

# Reference_based integration
MISC_NKT <- c(MISC_A1, MISC_A12, MISC_A2, MISC_A22, MISC_A3, MISC_A32, MISC_A4, MISC_A42, MISC_A52, MISC_A62, 
              MISC_HD1, MISC_HD2, MISC_HD3, MISC_HD4, MISC_HD5, MISC_HD6, MISC_HD7, MISC_HD8, MISC_HD9, MISC_HD10,
              MISC_HD11, MISC_HD12, MISC_HD13, MISC_CHD1, MISC_CHD2, MISC_CHD3, MISC_CHD4, MISC_CHD5, MISC_CHD6,
              MISC_P1, MISC_P2, MISC_P31, MISC_P32, MISC_P41, MISC_P42, MISC_P5, MISC_P6, MISC_P7)
for (i in 1:length(MISC_NKT)) {
  MISC_NKT[[i]] <- NormalizeData(MISC_NKT[[i]], verbose = TRUE)
  MISC_NKT[[i]] <- FindVariableFeatures(MISC_NKT[[i]], selection.method = "vst",
                                                   nfeatures = 2000, vecbose = TRUE)
}
MISC_NKT.anchors <- FindIntegrationAnchors(object.list = MISC_NKT, reference = c(13, 30), dims = 1:30)
MISC_NKT_integrated <- IntegrateData(anchorset = MISC_NKT.anchors, dims = 1:30)

DefaultAssay(MISC_NKT_integrated) <- "integrated"
MISC_NKT_integrated <- ScaleData(MISC_NKT_integrated)
MISC_NKT_integrated <- RunPCA(MISC_NKT_integrated, npcs = 30, veverbose = TRUE)

# check UMAP
DefaultAssay(MISC_NKT_integrated) <- "integrated"
MISC_NKT_integrated <- FindNeighbors(MISC_NKT_integrated, dims = 1:8)
MISC_NKT_integrated <- FindClusters(MISC_NKT_integrated, resolution = 0.9)
MISC_NKT_integrated <- RunUMAP(MISC_NKT_integrated, dims = 1:8)
p1 <- DimPlot(MISC_NKT_integrated, reduction = "umap", label = TRUE,  label.size = 3.5)
p1

# annotation
Idents(MISC_NKT_integrated) <- "integrated_snn_res.0.9"
DefaultAssay(MISC_NKT_integrated) <- "RNA"
MISC_NKT_integrated <- RenameIdents(MISC_NKT_integrated, '0' = "Naive CD4T cells_1", '1' = "Naive CD4T cells_2", '2' = "Naive CD8T cells", 
                                    '3' = "CD56dim CD38+ NK cells", '4' = "Memory CD4 T cells", '5' = "Terminal effector memory CD8 T cells", 
                                    '6' = "Effector memory CD8 T cells", '7' = "CCR6+ Memory CD4+ T cells", '8' = "Naive CD4 T cells_3", 
                                    '9' = "Vd2 gdT cells", '10' = "MAIT cells", '11' = "CD56bright NK cells", '12' = "CD56dim S1004A+ NK cells",
                                    '13' = "Regulatory T cells", '14' = "Central memory CD8  T cells", '15' = "CXCR3+ Memory CD4+ T cells",
                                    '16' = "T/NK/monocyte doublets", '17' = "Ki67+ T/NK cells")
p1 <- DimPlot(MISC_NKT_integrated, reduction = "umap", label = TRUE,  label.size = 3.5, repel = TRUE)
p1

# remove doublet_cluster
MISC_NKT <- subset(MISC_NKT_integrated, idents = "16", invert = T)

# saveRDS
saveRDS(MISC_NKT, "MISC_NKTclusters.rds")