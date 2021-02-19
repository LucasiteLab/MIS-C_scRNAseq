# Kenneth B. Hoehn
# 2/1/2021
# Combine BCR data from different sources into a single file
# annotated with metadata

library(alakazam)
library(dplyr)
library(stringr)
library(stringdist)

exclude = c(
"TP6A_BCR_VHT_cellranger","TP6B_BCR_VHT_cellranger", 
"TP7A_BCR_VHT_cellranger","TP7B_BCR_VHT_cellranger", 
"TS4A_BCR_VHT_cellranger","TS4B_BCR_VHT_cellranger", 
"TS5A_BCR_VHT_cellranger","TS5B_BCR_VHT_cellranger",
"34F-BCR_VHT_cellranger")

dirs = list.files("raw")
dirs = dirs[!dirs %in% exclude]

# read in metadata file
metadata = read.csv("processed/patient_info_table_final_edit.csv",
	stringsAsFactors=FALSE)

metadata[metadata$timepoint == "",]$timepoint = "A"

data = tibble()
for(dir in dirs){
	print(dir)
	h = readChangeoDb(paste0("raw/",
		dir,
		"/filtered_contig_heavy_productive-T.tsv"))
	l = readChangeoDb(paste0("raw/",
		dir,
		"/filtered_contig_light_productive-T.tsv"))

	comb = bind_rows(h,l)

	meta = filter(metadata,bcr_library_id == dir)

	# bind metadata columns to data
	if(nrow(meta) != 0){
		print(paste(dir,"matched metadata table"))
		comb$orig.barcode = comb$cell_id
		comb = bind_cols(comb,meta)
		comb$old.ident = comb$orig.ident
	}else{
		# figure out the sample information for CITE files
		print(paste(dir,"didn't match metadata table"))
		str = strsplit(dir,split="_")[[1]][1]
		citen = as.numeric(substr(str,5,5))
		citebarcodes = read.csv(paste0("cite/cite_seq",citen
			,"_kept_prop.csv"), stringsAsFactors=FALSE)
		print(paste(str,"barcode match",mean(comb$cell_id %in%
			citebarcodes$barcode)))
		# match sample IDs to cell barcodes in BCR data
		comb = filter(comb,cell_id %in% citebarcodes$barcode)
		m = match(comb$cell_id,citebarcodes$barcode)
		if(sum(comb$cell_id != citebarcodes$barcode[m]) != 0){
			stop("barcode match failure")
		}
		# match sample IDs to sample IDs in metadata file
		comb$sample = citebarcodes$sample[m]
		print(paste(paste(unique(comb$sample[!comb$sample %in%
			metadata$orig.ident]),
			collapse=" "),"not in metadata"))
		comb = filter(comb,sample %in% metadata$orig.ident)
		m = match(comb$sample,metadata$orig.ident)
		meta = metadata[m,]
		comb = select(comb,-sample)
		comb = bind_cols(comb,meta)
		comb$old.ident = paste0("CITE_",citen)
		comb$orig.barcode = paste("CITE",citen,comb$cell_id,sep="_")
	}
	data = bind_rows(data,comb)
}
table(data$sample_id,data$productive)
table(data$sample_id,data$old.ident)
data = filter(data,productive == TRUE)
writeChangeoDb(data,file="processed/alldata.tsv")

