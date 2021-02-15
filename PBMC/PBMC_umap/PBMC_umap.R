library(tidyverse)
library(Seurat)
library(ggplot2)
library(future)


## Import hashed COVID files


file_path <- "/gpfs/ysm/project/lucas/ar2374/COV/Data_COVID/"


# CITE samples

folders_cite <- c("CITE_1", "CITE_2", "CITE_3", "CITE_4", "CITE_5", "CITE_6") 

list_cite <- list()

for(i in 1:length(folders_cite)){
  tenx <- Read10X(data.dir= paste0(file_path, folders_cite[i], "/"))
  gex <- CreateSeuratObject(counts = tenx, min.cells = 5, project= folders_cite[i])
  list_cite[[i]] <- gex
}

names(list_cite) <- folders_cite


# COVID patient samples

folders_patient <- c("NS1A", "NS1B", "TS2A", "TS2B", "TS3A", "TS3B", "TP9B", "NS0A", "NS0B")

list_patient <- list()

for(i in 1:length(folders_patient)){
  tenX <- Read10X(data.dir= paste0(file_path, folders_patient[i], "/"))
  gex <- CreateSeuratObject(counts = tenX, min.cells = 5, project= folders_patient[i])
  list_patient[[i]] <- gex
}

names(list_patient) <- folders_patient


# Read in barcodes

prefix <- "Data_COVID/Dehashing/cite_seq" #have to be in right directory
num_string <- c( "1", "2", "3", "4", "5", "6" )
suffix <- "_kept_prop.csv"

hash_files <- list()
for(i in 1:length(num_string)) {
  file_name <- paste0(prefix, num_string[i], suffix)
  hash_cite <- read.csv(file_name)
  hash_cite[,4] <- rep(num_string[i], nrow(hash_cite))
  hash_files[[i]] <- hash_cite 
}

hash_cite1.2 <- rbind(hash_files[[1]], hash_files[[2]])
hash_cite3.4 <- rbind(hash_files[[3]], hash_files[[4]])
hash_cite5.6 <- rbind(hash_files[[5]], hash_files[[6]])

hash_files_list <- list()

hash_files_list[[1]] <- hash_cite1.2
hash_files_list[[2]] <- hash_cite1.2
hash_files_list[[3]] <- hash_cite3.4
hash_files_list[[4]] <- hash_cite3.4
hash_files_list[[5]] <- hash_cite5.6
hash_files_list[[6]] <- hash_cite5.6
hash_files_list[[7]] <- hash_cite5.6
hash_files_list[[8]] <- hash_cite3.4
hash_files_list[[9]] <- hash_cite3.4


# Find barcodes relevant to each patient and subset Seurat objects

list_subset <- list() 

cite1 <- c("CITE_1", "CITE_1", "CITE_3", "CITE_3", "CITE_5", "CITE_5", "CITE_5", "CITE_3", "CITE_3")

cite2 <- c("CITE_2", "CITE_2", "CITE_4", "CITE_4", "CITE_6", "CITE_6", "CITE_6", "CITE_4", "CITE_4")

hash <- c("3", "4", "5", "6", "3", "4", "8", "3", "4")

patient_lookup <- data.frame("patient" = folders_patient,
                             "cite1" = cite1,
                             "cite2" = cite2, 
                             "hash" = hash)

row.names(patient_lookup) <- patient_lookup[,1]

# make a list of subsetted Seurat objects and their correlates

for(i in 1:nrow(patient_lookup)) { 
  name <- patient_lookup[i,1]
  print(name)
  hash <- patient_lookup[i,4]
  pcite1 <- patient_lookup[i,2]
  pcite2 <- patient_lookup[i,3]
  hash_files <- hash_files_list[[i]]
  hash_subset <- hash_files %>% filter(hashing == paste0("#", hash) & sample == name) #filter df by patient's associated hash
  hash_vector <- unique(as.character(hash_subset[['barcode']]))  #set as character
  first_cite <- list_cite[[pcite1]][ , hash_vector]  #subset Seurat corresponding cite sample by relevant hash barcodes
  second_cite <- list_cite[[pcite2]][ , hash_vector]
  orig.ident1 <- as.character(first_cite@meta.data$orig.ident[1])
  orig.ident2 <- as.character(second_cite@meta.data$orig.ident[1])
  for(j in 1:length(hash_vector)){
    
    index <- grep(hash_vector[j], hash_files[,1], fixed = TRUE)
    if (length(index) > 1) {
      for (idx in 1:length(index)) {
        if (hash_files[index[idx], 3] != name) {
          barcode <- hash_files[index[idx], 1]
          orig.name <- paste0('CITE_', hash_files[index[idx], 4])
          if (orig.ident1 == orig.name) {
            first_cite <- subset(first_cite, cells = barcode, invert=TRUE)
          } else if (orig.ident2 == orig.name) {
            second_cite <- subset(second_cite, cells = barcode, invert=TRUE)
          } else {
            print('OH NOOOOOO!!!!')
          }
          
        }
      }
    }
  }
  
  first_cite$orig.ident <- NULL
  second_cite$orig.ident <- NULL
  first_cite$orig.ident <- rep(names(list_patient)[[i]], nrow(first_cite@meta.data)) #reassign orig.ident
  second_cite$orig.ident <- rep(names(list_patient)[[i]], nrow(second_cite@meta.data))
  tmp_list <- list()
  tmp_list[[1]] <- first_cite #store in first item of list, first index
  tmp_list[[2]] <- second_cite #store in second item of list, second index
  tmp_list[[3]] <- list_patient[[i]]
  names(tmp_list) <- c(pcite1, pcite2, names(list_patient)[[i]])
  list_subset[[i]] <- tmp_list
}

names(list_subset) <- patient_lookup$patient

# merge list items into individual elements by patient

list_by_patient <- list()


for(i in 1:length(list_subset)) {
  item <- list_subset[[i]]
  id_name <- names(list_subset)[i]
  print(id_name)
  seurat1 <- item[[1]]
  seurat2 <- item[[2]]
  seurat3 <- item[[3]]
  list_by_patient[[i]] <- merge(seurat1,  y = c(seurat2, seurat3), 
                                add.cell.ids = c(paste0(names(item)[1]), 
                                                 paste0(names(item)[2]), 
                                                 paste0(names(item)[3])),
                                project = id_name)
}

names(list_by_patient) <- patient_lookup$patient



## Import MISC and Healthy Donor files

file_path <- "/gpfs/ysm/project/lucas/ar2374/COV/Data_new/"


folders_misc_hd <- c("Y111-1", "Y112-1", "Y113-1", "HD_35F", "HD_32M", "HD_36M", 
                     "NC-13F", "Y117-1", "Y117-R", "Y124-1", "Y124-R", "Y125-1", 
                     "Y127-1", "Y129-1", "Y28-2", "Y28-4", "Y29-2", "Y54-4", 
                     "Y70-4",
                     "HA5876", "HA5877", "HA5894", "HA5952", "HA5957", "HA5953",
                     "C39","C32", "C27", "C33", "TP8B")


folders_misc_names <- c("Y111_1", "Y112_1", "Y113_1", "HD_35F", "HD_32M", "HD_36M", 
                        "NC_13F", "Y117_1", "Y117_R", "Y124_1", "Y124_R", "Y125_1", 
                        "Y127_1", "Y129_1", "Y28_2", "Y28_4", "Y29_2", "Y54_4", 
                        "Y70_4",
                        "HA5876", "HA5877", "HA5894", "HA5952", "HA5957", "HA5953",
                        "C39","C32", "C27", "C33", "TP8B")

# loop through directory  
for(i in 1:length(folders_misc_hd)){
  tenx <- Read10X(data.dir= paste0(file_path, folders_misc_hd[i], "/"))
  if(is.list(tenx)){
    tenx <- tenx[[1]]
  }
  gex <- CreateSeuratObject(counts = tenx, min.cells = 5, project= folders_misc_hd[i])
  list_by_patient[[9+i]] <- gex
  names(list_by_patient)[9+i] <- folders_misc_names[i]
}


saveRDS(list_by_patient, file = "list_by_patient/list_by_patient.rds")



## Data Processing 

list_by_patient <- readRDS("list_by_patient/list_by_patient.rds")


for (i in 1:length(list_by_patient)) {
  list_by_patient[[i]][["percent.mt"]] <- PercentageFeatureSet(list_by_patient[[i]], 
                                                               pattern = "^MT-")
  list_by_patient[[i]] <- subset(list_by_patient[[i]], 
                                 subset = nFeature_RNA > 200 & percent.mt < 10) 
  list_by_patient[[i]] <- NormalizeData(list_by_patient[[i]], 
                                        verbose = FALSE)
  list_by_patient[[i]] <- FindVariableFeatures(list_by_patient[[i]], 
                                               selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
}

list_by_patient[[11]] <- NULL

# be sure to set 4 cores, 500GB
plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2e10)
misc.anchors <- FindIntegrationAnchors(object.list = list_by_patient, # best umaps from reference based, dim30, regressed
                                       #                                       reduction = "rpca",
                                       reference = c(10,14), #modified to 14 for 112 exclusion
                                       dims = 1:30)



saveRDS(misc.anchors, file = "anchors_rds/misc_anchors_2ref_dim30_wo112.rds")

misc.int <- IntegrateData(anchorset = misc.anchors, dims = 1:30)


DefaultAssay(misc.int) <- "integrated"

saveRDS(misc.int, file = "integrated_rds/misc_int_2ref_dim30_wo112.rds")


# Cell cycle regression 


misc.int <- readRDS("integrated_rds/misc_int_2ref_dim30_wo112.rds") # or whatever integrated object

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

misc.int <- CellCycleScoring(misc.int, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#(Warning: The following features are not present in the object: DTL, PRIM1, MLF1IP, RPA2, CCNE2, UBR7, POLD3, MSH2, RAD51, CDC45, EXO1, BLM, CASP8AP2, POLA1, BRIP1, E2F8, not searching for symbol synonyms
#Warning: The following features are not present in the object: FAM64A, CKAP2L, GTSE1, HJURP, HN1, TTK, CDC25C, KIF2C, RANGAP1, DLGAP5, CDCA2, PSRC1, ANLN, LBR, CTCF, NEK2, G2E3, GAS2L3, CENPA, not searching for symbol synonyms)

plot1 <- RidgePlot(misc.int, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

ggsave(plot = plot1, file = "ridge_plot.png")

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2e10)
misc.int.regressed <- ScaleData(misc.int, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(misc.int))

saveRDS(misc.int.regressed, file = "/gpfs/ysm/project/lucas/ar2374/COV/integrated_rds/misc_int_regressed_dim30_wo112.rds")

# Rest of clustering


misc.int <- misc.int.regressed
misc.int <- RunPCA(misc.int, npcs = 35, verbose = FALSE) # should run on integrated
plot2 <- ElbowPlot(misc.int)
misc.cluster <- misc.int #save a copy with PCA done for UMAP and clustering later
ggsave("elbow_plots/pc_elbow_plot_2ref_dim30_regressed_wo112.pdf", plot = plot2)

misc.int <- RunUMAP(misc.int, reduction = "pca", dims = 1:30) #only run the first time

orig.ident <- as.vector(misc.int@meta.data$orig.ident)
Idents(misc.int) <- orig.ident

plot3 <- DimPlot(misc.int, reduction = "umap", group.by = "orig.ident") #only run the first time
ggsave("umap_batch_pbmc/umap_batch_30_covid_2ref_regressed_wo112.pdf", plot = plot3, height = 12, width = 14) #only run the first time


# Find neighbors and clusters (louvain) based on PCA, and then run umap 

misc.cluster <- FindNeighbors(misc.cluster, dims = 1:30) #not uMAP-ed data
misc.cluster <- FindClusters(misc.cluster, resolution = 1, random.seed= 10) #experiment with diff resolution
head(Idents(misc.cluster), 5)

misc.cluster <- RunUMAP(misc.cluster, dims = 1:30) #default is reduction = pca
plot4 <- DimPlot(misc.cluster, reduction = "umap", 
                 cols=DiscretePalette(36, palette='polychrome'),
                 pt.size =0.05, label = TRUE)
ggsave("umap_louvain_pbmc/umap_louvain_res1dim30_2ref_regressed_wo112.pdf", plot = plot4, height = 12, width = 14)

saveRDS(misc.cluster, file = "/gpfs/ysm/project/lucas/ar2374/COV/cluster_rds_pbmc/misc_cov_cluster_dim30_2ref_dim30_regressed_wo112.rds")

# Find Markers

library(future)
plan("multiprocess", workers = 8)
ref2_markers <- FindAllMarkers(misc.cluster, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

ref2_markers

write.csv(ref2_markers, file = "cluster_markers_pbmc/ref2_dim30_markers_regressed_wo112.csv")


# Cleaning the dataset

misc.cluster <- readRDS("cluster_rds_pbmc/priority_not_clean/misc_cov_cluster_dim30_2ref_dim30_regressed_wo112.rds") #full dataset integration with 112
list_by_patient <- readRDS("list_by_patient/list_by_patient.rds") # unmodified list_by_patient 


for (i in 1:length(list_by_patient)) {
  list_by_patient[[i]][["percent.mt"]] <- PercentageFeatureSet(list_by_patient[[i]], # initial filtering for re-clustering
                                                               pattern = "^MT-")
  list_by_patient[[i]] <- subset(list_by_patient[[i]], 
                                 subset = nFeature_RNA > 200 & percent.mt < 10) 
  list_by_patient[[i]] <- NormalizeData(list_by_patient[[i]], 
                                        verbose = FALSE)
  list_by_patient[[i]] <- FindVariableFeatures(list_by_patient[[i]], 
                                               selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
}

list_by_patient[[11]] <- NULL

query_vector <- rownames(misc.cluster@meta.data)
orig_ident <- misc.cluster@meta.data[,3]

query_vector_mod <- c()

for(i in 1:length(query_vector)){
  split_barcode <- strsplit(query_vector[i], "_")[[1]]
  if(length(split_barcode) == 1){
    query_vector_mod[i] <- split_barcode
  } else if(length(split_barcode) == 2){
    query_vector_mod[i] <- split_barcode[1]
  } else if(length(split_barcode) == 3){
    query_vector_mod[i] <- paste(split_barcode[1], split_barcode[2], sep = "_")
  } else if(length(split_barcode) == 4){
    query_vector_mod[i] <- paste(split_barcode[1], split_barcode[2], split_barcode[3], sep = "_")
  }
}



new_name_vec <- c()

orig_ident_names <- c("NS1A", "NS1B", "TS2A", "TS2B", "TS3A", "TS3B", "TP9B", "NS0A", "NS0B",
                      "Y111-1", "Y113-1", "HD_35F", "HD_32M", "HD_36M", 
                      "NC-13F", "Y117-1", "Y117-R", "Y124-1", "Y124-R", "Y125-1", 
                      "Y127-1", "Y129-1", "Y28-2", "Y28-4", "Y29-2", "Y54-4", 
                      "Y70-4",
                      "HA5876", "HA5877", "HA5894", "HA5952", "HA5957", "HA5953",
                      "C39","C32", "C27", "C33", "TP8B")


names(list_by_patient) <- orig_ident_names


for(i in 1:length(query_vector_mod)) {
  list_name <- orig_ident[i]
  new_name <- grep(query_vector_mod[i], rownames(list_by_patient[[list_name]]@meta.data), fixed = TRUE, value = TRUE)
  if(length(new_name) < 1) {
    print(query_vector_mod[i])
  } else {
    new_name_vec[i] <- new_name
  }
}

meta <- misc.cluster@meta.data
meta[,ncol(meta)+1] <- new_name_vec 

meta_remove <- meta %>% filter(seurat_clusters %in% c(15,18))

list_by_patient_clean <- list() # a new list for cleaned list_by_patient

for(i in 1:length(list_by_patient)) {
  to_remove <- meta_remove %>% filter(orig.ident == names(list_by_patient)[[i]])
  list_by_patient_clean[[i]] <- subset(list_by_patient[[i]], cells = to_remove[,11], invert = TRUE) # subset and populate clean list
  names(list_by_patient_clean)[[i]] <- names(list_by_patient)[[i]]
}

saveRDS(list_by_patient_clean, file = "list_by_patient/list_by_patient_clean_wo112.rds")


##

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2e10)
misc.anchors <- FindIntegrationAnchors(object.list = list_by_patient_clean, # best umaps from reference based, dim30, regressed
                                       #                                       reduction = "rpca",
                                       reference = c(10,14), #modified to 14 for 112 exclusion
                                       dims = 1:30)



saveRDS(misc.anchors, file = "anchors_rds/misc_anchors_2ref_dim30_regressed_wo112_clean.rds")

misc.int <- IntegrateData(anchorset = misc.anchors, dims = 1:30)


DefaultAssay(misc.int) <- "integrated"

saveRDS(misc.int, file = "integrated_rds/misc_int_2ref_dim30_regressed_wo112_clean.rds")


##

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

misc.int <- CellCycleScoring(misc.int, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#(Warning: The following features are not present in the object: DTL, PRIM1, MLF1IP, RPA2, CCNE2, UBR7, POLD3, MSH2, RAD51, CDC45, EXO1, BLM, CASP8AP2, POLA1, BRIP1, E2F8, not searching for symbol synonyms
#Warning: The following features are not present in the object: FAM64A, CKAP2L, GTSE1, HJURP, HN1, TTK, CDC25C, KIF2C, RANGAP1, DLGAP5, CDCA2, PSRC1, ANLN, LBR, CTCF, NEK2, G2E3, GAS2L3, CENPA, not searching for symbol synonyms)

plot1 <- RidgePlot(misc.int, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

ggsave(plot = plot1, file = "ridge_plot.png")

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2e10)
misc.int.regressed <- ScaleData(misc.int, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(misc.int))

saveRDS(misc.int.regressed, file = "/gpfs/ysm/project/lucas/ar2374/COV/integrated_rds/misc_int_regressed_2ref_dim30_wo112_clean.rds")


# Rest of clustering


misc.int <- misc.int.regressed
misc.int <- RunPCA(misc.int, npcs = 35, verbose = FALSE) # should run on integrated
plot2 <- ElbowPlot(misc.int)
misc.cluster <- misc.int #save a copy with PCA done for UMAP and clustering later
ggsave("elbow_plots/pc_elbow_plot_2ref_dim30_regressed_wo112_clean.pdf", plot = plot2)

misc.int <- RunUMAP(misc.int, reduction = "pca", dims = 1:30) #only run the first time

orig.ident <- as.vector(misc.int@meta.data$orig.ident)
Idents(misc.int) <- orig.ident

plot3 <- DimPlot(misc.int, reduction = "umap", group.by = "orig.ident") #only run the first time
ggsave("umap_batch_pbmc/umap_batch_30_covid_2ref_regressed_wo112_clean.pdf", plot = plot3, height = 12, width = 14) #only run the first time


# Find neighbors and clusters (louvain) based on PCA, and then run umap 

misc.cluster <- FindNeighbors(misc.cluster, dims = 1:30) #not uMAP-ed data
misc.cluster <- FindClusters(misc.cluster, resolution = 1, random.seed= 10) #experiment with diff resolution
head(Idents(misc.cluster), 5)

misc.cluster <- RunUMAP(misc.cluster, dims = 1:30) #default is reduction = pca
plot4 <- DimPlot(misc.cluster, reduction = "umap", 
                 cols=DiscretePalette(36, palette='polychrome'),
                 pt.size =0.05, label = TRUE)
ggsave("umap_louvain_pbmc/umap_louvain_res1dim30_2ref_regressed_wo112_clean.pdf", plot = plot4, height = 12, width = 14)

saveRDS(misc.cluster, file = "/gpfs/ysm/project/lucas/ar2374/COV/cluster_rds_pbmc/misc_cov_cluster_dim30_2ref_regressed_wo112_clean_res0.9.rds")


# Find markers

library(future)
plan("multiprocess", workers = 8)
ref2_markers <- FindAllMarkers(misc.cluster, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

ref2_markers

write.csv(ref2_markers, file = "cluster_markers_pbmc/markers_dim30_2ref_regressed_wo112_clean.csv")

  
# Sample Composition
  # This feeds into cell type distribution plots
cluster_prop <- table(Idents(misc.cluster), misc.cluster$orig.ident) #just the cells

write.csv(cluster_prop, file = "sample_full_integrated_pbmc_wo112.csv")


