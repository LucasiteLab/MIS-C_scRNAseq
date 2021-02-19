# Libraries
library(Seurat)
library(ggplot2)
library(AUCell)
library(dplyr)

# geneSets can take character vector or list if multiple signatures. 

# Extract expression matrix from Seurat object (monocyte)
misc.mono <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_updated_cond.rds")

# Subset Seurat; needs to be run individually due to size; AUCell calculated on the cell level
Idents(misc.mono) <- "condition_new"
cohort_idents <- c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

cohort_seurat <- list()
cohort_names <- list()

for(i in 1:length(cohort_idents)){
  cohort_seurat[[i]] <- subset(misc.mono, idents = cohort_idents[i])
  names(cohort_seurat)[i] <- cohort_idents[i]
  # cohort_names[[i]] <- colnames(cohort_seurat[[i]]) #only need the names if subsetting off of AUC object
  # names(cohort_names)[i] <- cohort_idents[i]
}

#cell_rank_list <- list()
cell_rank_list <- readRDS("AUCell/pbmc_rankings_by_cohort.rds")
cell_AUC_list <- list()

go_ifn <- read.csv("AUCell/AUC_go_ifn/go_response_type_I_ifn.csv", header = FALSE) # data frame consisting of gene names in GO: Response to type I IFN
go_ifn <- go_ifn[,1]

absent <- c() # remove gene names not also present in seurat objects
for(i in 1:length(go_ifn)){
  if(!any(rownames(cohort_seurat[[1]]@assays$RNA) == go_ifn[i])){
    absent <- c(absent, go_ifn[i])
  }
}

go_ifn_present <- setdiff(go_ifn, absent)

for(i in 1:length(cohort_seurat)){
  exprmat <- GetAssayData(cohort_seurat[[i]][['RNA']], slot = "counts")
  #  cell_rank_list[[i]] <- AUCell_buildRankings(exprmat)
  cell_AUC_list[[i]] <- AUCell_calcAUC(go_ifn_present, cell_rank_list[[i]], aucMaxRank = ceiling(0.1 * nrow(cell_rank_list[[i]]))) # set the top genes
}

saveRDS(cell_AUC_list, file = "AUCell/pbmc_AUC_by_cohort_go_ifn.rds")


# Plot histograms to explore gene expression data- do a large number of cells express lots of genes or the opposite? 

exprmat_list <- list()
for(i in 1:length(cohort_seurat)){ 
  exprmat <- GetAssayData(cohort_seurat[[i]][['RNA']], slot = "counts")
  exprmat_list[[i]] <- exprmat
  names(exprmat_list)[i] <- cohort_idents[i]
}

for(i in 1:length(exprmat_list)){ #exploring gene data
  countByCell_misc <- colSums(as.matrix(exprmat_list[[i]]), na.rm = TRUE)
  png(file = paste0("AUCell/", cohort_idents[i], "_pbmc_rank_histogram.png"))
  hist(countByCell_misc, xlim = c(0,8000), breaks = 100)
  dev.off()
}

# Transfer AUC values for each cohort to a data frame and then into Seurat obj

library(dplyr)
mod_seurat <- list()

for(i in 1:length(cell_AUC_list)){
  AUC_vec <- getAUC(cell_AUC_list[[i]])
  AUC_df <- as.data.frame(AUC_vec)
  AUC_df2 <- as.data.frame(t(as.matrix(AUC_df)))
  AUC_df2[,2] <- rownames(AUC_df2)
  metadata <- cohort_seurat[[i]]@meta.data
  metadata[,29] <- rownames(metadata)
  colnames(metadata)[29] <- "barcode_name"
  colnames(AUC_df2)[2] <- "barcode_name"
  metadata_AUC <- left_join(metadata, AUC_df2, by = "barcode_name")
  rownames(metadata_AUC) <- rownames(metadata)
  cohort_seurat[[i]]@meta.data <- metadata_AUC
  mod_seurat[[i]] <- cohort_seurat[[i]]
}

saveRDS(mod_seurat, file = "AUCell/pbmc_seurat_by_cohort_wAUC_go_ifn.rds")


# For monocytes

df_list <- list()
for(i in 1:length(mod_seurat)){
  Idents(mod_seurat[[i]]) <- "seurat_clusters"
  misc.mono <- subset(mod_seurat[[i]], idents = c(2,8,12)) #class monocytes, neutrophils, nonclass mono
  pbmc_meta <- misc.mono@meta.data
  df_auc <- data.frame('sample' = pbmc_meta$sample_id, 'score' = pbmc_meta$geneSet, 
                       'cluster' = pbmc_meta$annotation, 'condition' = pbmc_meta$condition_new)
  df_list[[i]] <- df_auc
  names(df_list)[i] <- cohort_idents[i]
}

summary_df_auc <- do.call(rbind, df_list)

write.csv(summary_df_auc, file = "AUCell/AUCell_myeloid_spec_csv/AUC_go_ifn_myeloid_specific.csv") #use to calculate pvalue 


# Explore AUC thresholds

summary_df_auc <- read.csv("AUCell/AUC_bac1_myeloid.csv") # or pbmc

myeloid_hist <- ggplot(summary_df_auc, aes(x = score)) + #histogram
  geom_histogram(binwidth = 0.01, color = "black", fill = "lightblue") +
  # geom_vline(xintercept = 0.165, linetype = "longdash")+
  ylab("Frequency")+
  xlab("AUC")+
  theme_classic(base_size = 16)
ggsave(myeloid_hist, file = "AUCell/myeloid_spec_hist_AUC_go_ifn.png", height = 5, width = 5)


# Box plot
library(ggplot2)
library(dplyr)

cyto <- read.csv("Sheets/AUC_go_ifn_monocyte_specific.csv")

names <- unique(cyto$sample)
means <- data.frame("name" = rep(NA,38), "value" = rep(NA,38))

for(i in 1:length(names)){
  cyto_tmp <- cyto %>% filter(sample == names[i])
  means[i,1] <- names[i]
  means[i,2] <- mean(cyto_tmp[,3])
}

means[means$name %in% c("A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5", "A.HD6",
                        "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11",
                        "A.HD12", "A.HD13"), 'condition'] <- 'A.HD'

means[means$name %in% c("A.COV1.1", "A.COV2.1", "A.COV3.1","A.COV4.1"), 'condition'] <- 'COVID19-A'

means[means$name %in% c("A.COV1.2", "A.COV2.2", "A.COV3.2", "A.COV4.2",
                        "A.COV5.2", "A.COV6.2"), 'condition'] <- 'COVID19-B'

means[means$name %in% c("P1.1", "P2.1", "P3.1", "P4.1",
                        "P5.1", "P6.1", "P7.1"), 'condition'] <- 'MIS-C'

means[means$name %in% c("C.HD1", "C.HD2", "C.HD3", 
                        "C.HD4", "C.HD5", "C.HD6"), 'condition'] <- 'C.HD'

means[means$name %in% c("P3.2", "P4.2"), 'condition'] <- 'MIS-C-R'

means_severe <- means %>% filter(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1"))


level_order <-  c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

means$condition <- factor(means$condition, level = level_order)
means_severe$condition <- factor(means_severe$condition, level = level_order)

misc_tmp <- means %>% filter(condition == "MIS-C")
chd_tmp <- means %>% filter(condition == "C.HD")

covid19a_tmp <- means %>% filter(condition == "COVID19-A")
covid19b_tmp <- means %>% filter(condition == "COVID19-B")
ahd_tmp <- means %>% filter(condition == "A.HD")

# Example

w.test <- wilcox.test(x = chd_tmp$value, y = misc_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval <- w.test$p.value


means2 <- means %>% filter(!(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))

cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

plot1 <- ggplot(means, aes(x = condition, y = value)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) + 
  geom_jitter(data=means_severe, colour ="#c94040",  size = 0.5, width = 0.1)+
  geom_jitter(data = means2, aes(colour = factor(condition, level = level_order)), size = 0.5, width = 0.1) +
  scale_color_manual(values = cols) +
  ggtitle("GO: Response to type I IFN") +
  xlab("") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylab("AUC") 

plot1

ggsave(plot1, file = "Fig2_viral_bac_fig/GO_response_IFN_AUC_sig_myeloid_specific.pdf", height = 1.5, width =2)




