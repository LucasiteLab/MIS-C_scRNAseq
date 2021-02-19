library(Seurat)
library(ggplot2)
library(dplyr)

# or any other cluster file with some samples containing ADT layer
misc.cluster <- readRDS("cluster_rds_pbmc/dim20/misc_cov_cluster_2ref.rds") 

misc.idents <- misc.cluster
Idents(misc.idents) <- "orig.ident"

## Make a smaller Seurat object with just the relevant cells from the six samples

# Make reference vector of barcodes

file_path <- "/gpfs/ysm/project/lucas/ar2374/COV/Data_new/"


folders_misc_hd <- c("Y111-1", "Y113-1", 
                     "HD_32M", "HD_35F", "HD_36M")

misc_hd_names <- c("Y111_1", "Y113_1",
                   "HD_32M", "HD_35F", "HD_36M")

hiro.names <- c("MISC_P1", "MISC_P2", "MISC_HD1", "MISC_HD2", "MISC_HD3")

barcodes_adt <- list()

for(i in 1:length(folders_misc_hd)) {
  barcodes_adt[[i]] <- WhichCells(misc.idents, idents = c(hiro.names[i]))
}

list_cite <- list()

# Read in ADT data

for(i in 1:length(folders_misc_hd)){
  tenx <- Read10X(data.dir= paste0(file_path, folders_misc_hd[i], "/", "CITE"))
  
  rownames(x = tenx[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
                                                   x = rownames(x = tenx[["Antibody Capture"]]))
  adt <- CreateSeuratObject(counts = tenx[["Antibody Capture"]], min.cells = 5, project = folders_misc_hd[i])
  list_cite[[i]] <- adt
  names(list_cite)[i] <- misc_hd_names[i]
}


# Create a vector of ADT-matching barcodes from larger dataset

barcodes_vector <- c()

for(i in 1:length(list_cite)){
  barcodes_mod <- c()
  col_names <- colnames(list_cite[[i]])
  
  for(j in 1:length(col_names)){
    anchor_barcode <- strsplit(col_names[j], '-', fixed = TRUE)[[1]][1]
    new_name <- grep(anchor_barcode, barcodes_adt[[i]], fixed = TRUE, value = TRUE)
    #print(new_name)
    if(length(new_name) > 0 ){ #used to be 0
      print(new_name)
      barcodes_mod[j] <- new_name #some barcodes may not be present in bigger data bc filtering
    } # } else {
    #   barcodes_mod[j] <- anchor_barcode #the barcodes that are kept the same will be removed later
    # }
  }
  barcodes_vector <- c(barcodes_vector, barcodes_mod[!is.na(barcodes_mod)])
  #list_cite[[i]] <- RenameCells(list_cite[[i]], new.names = barcodes_mod)
}

misc.six <- subset(misc.cluster, cells = barcodes_vector) # keep only these cells 


## Modify names of ADT objects themselves to add as an assay to subsetted seurat

# Read in as before

list_cite <- list()


for(i in 1:length(folders_misc_hd)){
  tenx <- Read10X(data.dir= paste0(file_path, folders_misc_hd[i], "/", "CITE"))
  
  rownames(x = tenx[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
                                                   x = rownames(x = tenx[["Antibody Capture"]]))
  adt <- CreateSeuratObject(counts = tenx[["Antibody Capture"]], min.cells = 5, project = folders_misc_hd[i])
  list_cite[[i]] <- adt
  names(list_cite)[i] <- misc_hd_names[i]
}


to_remove <- c()
list_cite_clean <- list()


for(i in 1:length(list_cite)){
  barcodes_mod <- c()
  col_names <- colnames(list_cite[[i]])
  
  for(j in 1:length(col_names)){
    anchor_barcode <- strsplit(col_names[j], '-', fixed = TRUE)[[1]][1]
    new_name <- grep(anchor_barcode, barcodes_adt[[i]], fixed = TRUE, value = TRUE)
    print(new_name)
    if(length(new_name) > 0 ){
      barcodes_mod[j] <- new_name #some barcodes may not be present in bigger data bc filtering
    } else {
      barcodes_mod[j] <- anchor_barcode #the barcodes that are kept the same will be removed later
      to_remove[j] <- anchor_barcode
    }
  }
  #barcodes_vector <- c(barcodes_vector, barcodes_mod[!is.na(barcodes_mod)])
  list_cite[[i]] <- RenameCells(list_cite[[i]], new.names = barcodes_mod)
  list_cite_clean[[i]] <- subset(list_cite[[i]], cells = to_remove, invert= TRUE)
}

for(i in 1:length(list_cite)){
  print(length(colnames(list_cite[[i]])))
}

# Now i have a clean list of ADT Seurat objects. I should be able to add these back to the bigger data

# Merge

merge.adt <- merge(x = list_cite_clean[[1]],
                   y = list_cite_clean[2:length(list_cite_clean)],
                   add.cell.ids = names(list_cite_clean), 
                   project = "Merge.Adt")

merge.adt[["percent.mt"]] <- PercentageFeatureSet(merge.adt, pattern = "^MT-")


# Add to data 

misc.six[["ADT"]] <- CreateAssayObject(merge.adt@assays$RNA@data)


misc.six <- NormalizeData(misc.six, assay = "ADT")


plot9 <- FeaturePlot(object = misc.six, 
                     features = "0171-anti-humanmouserat-CD278",
                     min.cutoff = "q05", 
                     max.cutoff = "q95",
                     pt.size = 1,
 #                    cols = c("lightgrey", "#ED6925"), 
                     order = TRUE) + theme_classic(base_size = 36)+ 
                    scale_colour_gradient2(low = "#000004", mid = "#BB3754", high = "#FCFFA4",
                         limits=c(1.5, 4.5), midpoint = 3) +
                    theme(axis.text = element_text(size = 25),
                    legend.text = element_text(size = 20),
                    plot.title = element_text(size = 28, face= "italic", hjust = 0.5)) +
                    ggtitle("CD278-ADT") 

ggsave("Figure_tcell/CD278_adt_inferno.pdf", plot = plot9, height = 10, width = 10)


