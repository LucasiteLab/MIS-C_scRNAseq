library(Seurat)
library(SeuratData)
library(Connectome)
library(ggplot2)
library(cowplot)

## Connectome with downsampling and integration 

# If doing fine annotations, use this step
# misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0.rds")
# misc.myeloid <- readRDS("shared/Myeloid/misc_integrated_myeloid_updated_cond.rds")
# misc.tcell <- readRDS("shared/t_cell/misc_integrated_tcell_updated_cond.rds")
# misc.bcell <- readRDS("shared/b_cell/misc_integrated_object_2.0_bcell_final.rds")
# 
# pbmc_meta <- misc.cluster@meta.data
# myeloid_meta <- misc.myeloid@meta.data
# tcell_meta <- misc.tcell@meta.data
# bcell_meta <- misc.bcell@meta.data
# 
# pbmc_meta$cell_barcode <- rownames(pbmc_meta)
# myeloid_meta$cell_barcode <- rownames(myeloid_meta)
# tcell_meta$cell_barcode <- rownames(tcell_meta)
# bcell_meta$cell_barcode <- rownames(bcell_meta)
# 
# tcell_info <- data.frame("cell_barcode" = tcell_meta$cell_barcode, 
#                          "fine_annotation" = tcell_meta$subcluster_annotation1)
# 
# myeloid_info <- data.frame("cell_barcode" = myeloid_meta$cell_barcode, 
#                          "fine_annotation" = myeloid_meta$FinalAnnotation)
# 
# bcell_info <- data.frame("cell_barcode" = bcell_meta$cell_barcode, 
#                            "fine_annotation" = bcell_meta$annotation_bcell)
# 
# merged_info <- rbind(tcell_info, myeloid_info, bcell_info)
# 
# pbmc_tcell <- left_join(pbmc_meta, merged_info, by = "cell_barcode")
# pbmc_updated <- pbmc_tcell
# rownames(pbmc_updated) <- rownames(pbmc_meta)
# 
# misc.cluster@meta.data <- pbmc_updated
# 
# saveRDS(misc.cluster, "connectome/pbmc_fine_annotation.rds")

# Start analysis - here did not use fine annotation because of limited cell numbers per parcel

misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0.rds")
Idents(misc.cluster) <- "condition"
misc.misc <- subset(misc.cluster, idents = c("MIS-C"))
misc.chd <- subset(misc.cluster, idents = c("C.HD"))

# Check the cluster sums: start with PBMC object and MISC

cluster_vec <- unique(misc.misc$annotation)
cohort_vec <- unique(misc.misc$condition)
sample_vec <- unique(misc.misc$sample_id)

subject_list <- list() # list of subject seurat object with names
cluster_int <- list() # list of vectors for cluster intersection
Idents(misc.misc) <- "sample_id"


for(i in 1:length(sample_vec)){
  subject_list[[i]] <- subset(misc.misc, idents = sample_vec[i])
  names(subject_list)[i] <- sample_vec[i]
  cluster_int[[i]] <- unique(subject_list[[i]]$annotation)
  names(cluster_int)[i] <- sample_vec[i]
}

cluster_vec_new <- Reduce(intersect, cluster_int)

container_misc <- matrix(nrow=length(cluster_vec_new), ncol=length(sample_vec))
container_chd <- matrix(nrow=length(cluster_vec_new), ncol=length(sample_vec))

container_misc <- as.data.frame(container_misc)
container_chd <- as.data.frame(container_chd)

rownames(container_misc) <- cluster_vec_new
colnames(container_misc) <- sample_vec



for(i in 1:length(cluster_vec_new)){ # Naive B-cells 
  for(j in 1:length(subject_list)){ # P1.1 
    anntable <- table(subject_list[[j]]$annotation)
    print(cluster_vec_new[i])
    num <- anntable[[cluster_vec_new[i]]]
    container_misc[i,j] <- num
  }
}


# Check the cluster sums: start with PBMC object and CHD

cluster_vec <- unique(misc.chd$annotation)
cohort_vec <- unique(misc.chd$condition)
sample_vec <- unique(misc.chd$sample_id)

subject_list <- list() # list of subject seurat object with names
cluster_int <- list() # list of vectors for cluster intersection
Idents(misc.chd) <- "sample_id"


for(i in 1:length(sample_vec)){
  subject_list[[i]] <- subset(misc.chd, idents = sample_vec[i])
  names(subject_list)[i] <- sample_vec[i]
  cluster_int[[i]] <- unique(subject_list[[i]]$annotation)
  names(cluster_int)[i] <- sample_vec[i]
}

cluster_vec_new <- Reduce(intersect, cluster_int)

container_chd <- matrix(nrow=length(cluster_vec_new), ncol=length(sample_vec))
container_chd <- as.data.frame(container_chd)

rownames(container_chd) <- cluster_vec_new
colnames(container_chd) <- sample_vec


for(i in 1:length(cluster_vec_new)){ # Naive B-cells 
  for(j in 1:length(subject_list)){ # P1.1 
    anntable <- table(subject_list[[j]]$annotation)
    print(cluster_vec_new[i])
    num <- anntable[[cluster_vec_new[i]]]
    container_chd[i,j] <- num
  }
}

# Check cohort sums and remove those less than 75

container_misc$sum_misc <- rep(NA, nrow(container_misc))
for(i in 1:nrow(container_misc)){
  ncols <- ncol(container_misc)
  container_misc[i,ncols] <- sum(container_misc[i,1:ncols-1])
}

container_chd$sum_chd <- rep(NA, nrow(container_chd))
for(i in 1:nrow(container_chd)){
  ncols <- ncol(container_chd)
  container_chd[i,ncols] <- sum(container_chd[i,1:ncols-1])
}

container_misc_filtered <- container_misc %>% filter(sum_misc > 75) #none removed
container_chd_filtered <- container_chd %>% filter(sum_chd > 75)


merged_numbers <- merge(container_misc_filtered , container_chd_filtered, by="row.names", 
                        all.x=TRUE, all.y = TRUE)
merged_numbers <- na.omit(merged_numbers)

#merged_numbers$sum_chd <- NULL
#merged_numbers$sum_misc <- NULL

rownames(merged_numbers) <- merged_numbers[,1]
merged_numbers[,1] <- NULL
merged_numbers$min_cohort <- rep(NA, nrow(merged_numbers))
for(i in 1:nrow(merged_numbers)){
  ncols <- ncol(merged_numbers)
  merged_numbers[i,ncols] <- min(merged_numbers[i,8], merged_numbers[i,15])
}

merged_numbers$annotation <- rownames(merged_numbers)
#doublets <- c("T/NK/monocyte_doublets")
doublets <- c("T-NK-B cell doublets", "Platelet-T cell doublets")
merged_no_db <- merged_numbers %>% filter(!(annotation %in% doublets))
merged_no_db$annotation <- NULL


# Downsampling

condition_all <- c("MIS-C", "C.HD")
subject_all <- list()

Idents(misc.cluster) <- "condition"
for(i in 1:length(condition_all)){
  subject_all[[i]] <- subset(misc.cluster, idents = condition_all[i])
  names(subject_all)[i] <- condition_all[i]
}

cluster_names <- rownames(merged_no_db)
min_cells <- merged_no_db$min_cohort

lookup <- data.frame("cluster" = cluster_names, "min" = min_cells)


downsampled <- list()

for(i in 1:length((subject_all))){
  Idents(subject_all[[i]]) <- "annotation"
}

for(i in 1:length(subject_all)){ # will run through every subject
  barcodes <- c()  # empty barcode for subsetting
  
  for(j in 1:nrow(lookup)){ # will run through every cluster
    barcodes_tmp <- c() # empty barcode for sampling
    subject_cells <- WhichCells(subject_all[[i]], idents = lookup[j,1])
    
    if(lookup[j,2] > 30) {
      barcodes_tmp <- sample(subject_cells, size = lookup[j,2], replace = FALSE)
      
    } else {
      print(lookup[j,1])
      barcodes_tmp <- subject_cells
    }
    barcodes <- c(barcodes, barcodes_tmp)
  }
  downsampled[[i]] <- subject_all[[i]][ , barcodes] #subset with accumulated barcodes
  names(downsampled)[[i]] <- names(subject_all)[[i]]
}

saveRDS(downsampled, "connectome/list_dn_by_cohort_broad.rds")



## Use immunological genes for connectome 

misc.dn <- merge(x = downsampled[[1]],
                 y = downsampled[2:length(downsampled)],
                 add.cell.ids = names(downsampled), 
                 project = "Downsampled")

Idents(misc.dn) <- "annotation"
DefaultAssay(misc.dn) <- "RNA"

immune_genes <- read.csv("connectome/connectome_genes_immune.csv") 
immune_genes[,4] <- NULL
immune.conn.genes <- union(immune_genes$Ligand, immune_genes$Receptor)
genes <- immune.conn.genes[immune.conn.genes %in% rownames(misc.dn)]

connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes_basic <- connectome.genes[connectome.genes %in% rownames(misc.dn)]

# Merge downsampled

for(i in 1:length((downsampled))){
  Idents(downsampled[[i]]) <- "annotation"
  DefaultAssay(downsampled[[i]]) <- "RNA"
}

misc.list <- downsampled

misc.con.list <- list()
for (i in 1:length(misc.list)){
  misc.list[[i]] <- ScaleData(misc.list[[i]],features = rownames(misc.list[[i]]))
  misc.con.list[[i]] <- CreateConnectome(misc.list[[i]],
                                         species = 'human',p.values = F, 
                                         min.cells.per.ident = 75, custom.list = immune_genes)
}

names(misc.con.list) <- names(misc.list)

diff <- DifferentialConnectome(misc.con.list[[2]],misc.con.list[[1]])



celltypes <- as.character(unique(Idents(misc.dn)))
celltypes.chd <- paste(celltypes, 'C.HD', sep = '_')
celltypes.misc <- paste(celltypes, 'MIS-C', sep = '_')
misc.dn$celltype.condition <- paste(Idents(misc.dn), misc.dn$condition, sep = "_")
misc.dn$celltype <- Idents(misc.dn)
Idents(misc.dn) <- "celltype.condition"

#
diff.p <- data.frame()
for (i in 1:length(celltypes)){
  temp <- FindMarkers(misc.dn,
                      ident.1 = celltypes.misc[i],
                      ident.2 = celltypes.chd[i],
                      verbose = FALSE,
                      features = genes,
                      min.pct = 0,
                      logfc.threshold = 0.1)  
  temp2 <- subset(temp, p_val_adj < 0.05)
  if (nrow(temp2)>0){
    temp3 <- data.frame(genes = rownames(temp2),cells = celltypes[i])
    diff.p <- rbind(diff.p, temp3)
  }
}
diff.p$cell.gene <- paste(diff.p$cells, diff.p$genes,sep = '.')

diff$source.ligand <- paste(diff$source,diff$ligand,sep = '.')
diff$target.receptor <- paste(diff$target,diff$receptor,sep = '.')
diff.sub <- subset(diff, source.ligand %in% diff.p$cell.gene & target.receptor %in% diff.p$cell.gene)

diff.up <- diff.sub
level.order <- c("CD4 naive I", "CD4 naive II", "CD8 naive", 
                 "Regulatory T cells", "CD8 memory", "CD4 and CD8 activated memory", 
                 "CD8 effector", "CD56dim CD16bright NK", "CD56bright CD16dim NK", 
                 "MAIT-NKT cells", "gdT cells", "Neutrophils", "Classical monocytes", 
                 "Non classical monocytes", "conventional DC", "Naive B", "Memory B", 
                 "Plasma cells")

diff.up$source <- factor(diff.up$source, levels = level.order)

write.csv(diff.sub, "connectome/diff_connectome_LFC0.1_miscVctrl.csv")

pdf("connectome/dn_cohort_broad_DSP.pdf", height = 10, width = 30)
DifferentialScoringPlot(diff.sub,min.score = 1,min.pct = 0.1,infinity.to.max = T)
dev.off()

diff.up.up <- subset(diff.up,ligand.norm.lfc > 1 & recept.norm.lfc > 1 )
pdf("connectome/dn_cohort_broad_LgUP_RcUP_ordered_LR1.pdf", height = 17, width = 18)
CircosDiff(diff.up.up,min.score = 2,min.pct = 0.1,lab.cex = 1)
dev.off()





