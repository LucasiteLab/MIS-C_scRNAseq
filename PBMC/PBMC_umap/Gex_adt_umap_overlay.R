# Libraries
library(Seurat)
library(ggplot2)
library(AUCell)
library(dplyr)


misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_updated_cond.rds")


## Defining different populations

# NCAM1

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_NCAM1",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("NCAM1")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_NCAM1.pdf", plot = plot9, height = 5, width = 5, units = "in")


# CD3D

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_CD3D",
                     min.cutoff = 1,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2.7) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD3D")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_CD3D.pdf", plot = plot9, height = 5, width = 5, units = "in")

# S100A8

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_S100A8",
                     min.cutoff = 1,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("S100A8")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_S100A8.pdf", plot = plot9, height = 5, width = 5, units = "in")


# PBPP

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_PPBP",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 4.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("PPBP")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_PPBP.pdf", plot = plot9, height = 5, width = 5, units = "in")


# PBPP

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_CD8A",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD8A")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_CD8A.pdf", plot = plot9, height = 5, width = 5, units = "in")


# MKI67

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_MKI67",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("MKI67")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_MKI67.pdf", plot = plot9, height = 5, width = 5, units = "in")


# CD14

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_CD14",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2.25) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD14")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_CD14.pdf", plot = plot9, height = 5, width = 5, units = "in")

# MS4A1

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_MS4A1",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("MS4A1")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_MS4A1.pdf", plot = plot9, height = 5, width = 5, units = "in")

# FCGR3A

plot9 <- FeaturePlot(object = misc.cluster, 
                     features = "rna_FCGR3A",
                     min.cutoff = 0,
                     pt.size = 0.75,
                     order = TRUE) + 
  scale_color_gradient2(low = "#cccccc", mid = "#cc6868", high = "#d60606", midpoint = 2.25) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("FCGR3A")

ggsave("RC_umaps/gex_pbmc/gex_pbmc_FCGR3A.pdf", plot = plot9, height = 5, width = 5, units = "in")


# PBMC cite seq

misc.adt <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_adt.rds")

misc.adt <- NormalizeData(misc.adt, assay = "ADT")

# CD66b

inferno_mod <- inferno(20)[3:18]

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0166-anti-human-CD66b",
#                     min.cutoff = "q05", 
                     min.cutoff = 0,
                     max.cutoff = "q99",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD66b")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc_new/adt_pbmc_CD66b.pdf", plot = plot9, height = 5, width = 5, units = "in")


#CD45ra

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0063-anti-human-CD45RA",
#                     min.cutoff = "q05", 
                     min.cutoff = "0", 
                     max.cutoff = "q99",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3.25) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD45RA")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc_new/adt_pbmc_CD45RA.pdf", plot = plot9, height = 5, width = 5, units = "in")


#CD45RO

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0087-anti-human-CD45RO",
                     min.cutoff = 0, 
#                     max.cutoff = "q99",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD45RO")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc_new/adt_pbmc_CD45RO.pdf", plot = plot9, height = 5, width = 5, units = "in")


#CD71

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0394-anti-human-CD71",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2.25) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD71")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc/adt_pbmc_CD71.pdf", plot = plot9, height = 5, width = 5, units = "in")


#HLADR

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0159-anti-human-HLADR",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 4) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("HLADR")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc/adt_pbmc_HLADR.pdf", plot = plot9, height = 5, width = 5, units = "in")


#HLADR

plot9 <- FeaturePlot(object = misc.adt, 
                     features = "0168-anti-human-CD57-R",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.3,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 4) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
  ggtitle("CD57")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.3
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_pbmc/adt_pbmc_CD57.pdf", plot = plot9, height = 5, width = 5, units = "in")


## T cell

misc.adt <- readRDS("shared/t_cell/misc_integrated_tcell_adt.rds")

misc.adt <- NormalizeData(misc.adt, assay = "ADT")

# CD4

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0072-anti-human-CD4",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) +
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD4")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_tcell_new/adt_tcell_CD4.pdf", plot = plot9, height = 5, width = 5, units = "in")


# CD38

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0389-anti-human-CD38",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD38")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_tcell_new/adt_tcell_CD38.pdf", plot = plot9, height = 5, width = 5, units = "in")


# HLADR

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0159-anti-human-HLADR",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3.2) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("HLA-DR")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_tcell_new/adt_tcell_HLADR.pdf", plot = plot9, height = 5, width = 5, units = "in")



# CD223

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0152-anti-human-CD223",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD223")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_tcell_new/adt_tcell_CD223.pdf", plot = plot9, height = 5, width = 5, units = "in")


# CD45RO

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0087-anti-human-CD45RO",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2.75) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD45RO")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_tcell_new/adt_tcell_CD45RO.pdf", plot = plot9, height = 5, width = 5, units = "in")


## B cell supplement

# CD278

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0171-anti-humanmouserat-CD278",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("ICOS")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_s5/adt_tcell_CD278_new.pdf", plot = plot9, height = 5, width = 5, units = "in")


# CD279

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0088-anti-human-CD279",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 3.25) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("PD-1")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_s5/adt_tcell_CD279_new.pdf", plot = plot9, height = 5, width = 5, units = "in")

# CD71

plot9 <- FeaturePlot(object = misc.adt, 
                     features =  "0394-anti-human-CD71",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD71")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_s5/adt_tcell_CD71_new.pdf", plot = plot9, height = 5, width = 5, units = "in")



# CD95

plot9 <- FeaturePlot(object = misc.adt, 
                     features = "0156-anti-human-CD95",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 0.5,
                     order = TRUE) + 
  scale_color_gradient2(low = inferno_mod[2], mid = "#C43C4EFF", high = "#F1ED6FFF", midpoint = 2.5) +
  theme_classic()+
  theme(axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 28, hjust = 0.5)) +
  ggtitle("CD95")

q <- ggplot_build(plot9)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.5
q$data[[1]]$stroke <- 0
#q$data[[2]]$alpha <- 0.5
plot9 <- ggplot_gtable(q)

ggsave("RC_umaps/adt_s5/adt_tcell_CD95_new.pdf", plot = plot9, height = 5, width = 5, units = "in")


