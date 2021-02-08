# Kenneth B. Hoehn
# 2/1/2021
# Match cells to barcodes and annotations
# Cluster sequences into clonal clusters
# Reconstruct germline sequences

library(alakazam)
library(shazam)
library(dplyr)
library(tidyr)

data = readChangeoDb(file="processed/alldata.tsv")
odata = data
data[data$sex == "F",]$sex = "Female"
data[data$sex == "M",]$sex = "Male"

bcells = tibble(read.csv("processed/misc_integrated_bcell_barcodes_2.1.csv"
    ,stringsAsFactors=FALSE))

# sanity check annotations and shared sample_ids
print(unique(data$sample_id[data$subject_id %in% bcells$subject_id]))
print(unique(data$sample_id[!data$subject_id %in% bcells$subject_id]))
print(table(bcells$annotation_bcell))

# make column of unique barcode:sample combinations
bcells$sample_idcode = paste0(bcells$raw.barcode,"_",bcells$old.ident)
data$sample_idcode = paste0(data$cell_id,"_",data$old.ident)
max(table(bcells$sample_idcode))
#max(table(bcells$sample_idcode_raw))

data = filter(data,sample_id %in% bcells$sample_id)
data$id_file = paste0(data$sequence_id,data$bcr_library_id)

#check that we aren't dropping any samples
print(mean(data$sample_idcode %in% bcells$sample_idcode))
print(mean(bcells$sample_idcode %in% data$sample_idcode))
print(table(data[data$sample_idcode %in%
 bcells$sample_idcode,]$sample_id))
print(table(data[!data$sample_idcode %in%
 bcells$sample_idcode,]$sample_id))

#filter to only cells with barcodes in GEX data
fdata = filter(data,sample_idcode %in%
    bcells$sample_idcode)

#check that we didn't lose all the CITE cells
print(length(fdata$sample_idcode[grepl("CITE",
  fdata$sample_idcode)]))

f = match(fdata$sample_idcode,bcells$sample_idcode)
print(mean(fdata$cell_id == bcells[f,]$raw.barcode))

# add annotations
fdata$cell_type = bcells[f,]$annotation
fdata$bcell_type = bcells[f,]$annotation_bcell

#subset to heavy and light chains
h = filter(fdata,locus=="IGH")
l = filter(fdata,locus!="IGH")

# Reconstruct germline sequences 
hrefs = c("~/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta",
	"~/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta",
	"~/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta")
exec = "~/.local/bin/CreateGermlines.py"

writeChangeoDb(h,"processed/heavy_chain.tsv")

#path.expand(exec)
r <- paste(path.expand(hrefs), collapse = " ")
command <- paste("CreateGermlines.py -d", "processed/heavy_chain.tsv", "-r", r, "--format airr",
	"--outname heavy -g dmask")
system(command)

h_germline = readChangeoDb("processed/heavy_germ-pass.tsv")

# Make distance to nearest sequence neighbor plots to identify
# optimal clonal distance threshold
dist_cross <- distToNearest(h_germline, sequenceColumn="junction", 
                        vCallColumn="v_call", jCallColumn="j_call",
                        model="ham", normalize="len", nproc=1,
                        cross="sample_id")

pdf("intermediates/crossDistance.pdf",height=20,width=8)
ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
             aes(x=cross_dist_nearest)) + 
    theme_bw() + 
    xlab("Cross-sample_id Hamming distance") + 
    ylab("Count") +
    geom_histogram(color="white", binwidth=0.02) +
    geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
    facet_grid(sample_id ~ ., scales="free_y")
dev.off()

plots = list()
pdf("intermediates/dist_to_nearest.pdf",width=6,height=6)
for(sample_id in unique(dist_cross$sample_id)){
   print(sample_id)
   temp = filter(dist_cross,!!sample_id==sample_id)
   dist_ham <- distToNearest(temp, sequenceColumn="junction", 
                       vCallColumn="v_call", jCallColumn="j_call",
                       model="ham", normalize="len", nproc=2)
   output <- findThreshold(dist_ham$dist_nearest, method="density")
   threshold <- output@threshold
   g = ggplot(subset(dist_ham, !is.na(dist_nearest)),aes(x=dist_nearest,
   ,y = ..density..)) + 
     theme_bw() + 
     xlab("Hamming distance") + 
     ylab("Count") +
     scale_x_continuous(breaks=seq(0, 1, 0.1)) +
     geom_histogram(color="white", binwidth=0.02) +
     ggtitle(paste(sample_id,threshold))+
     geom_histogram(
       aes(x=cross_dist_nearest,y = -..density..),
       color="white", binwidth=0.02,fill="black")+
     xlim(0,max(filter(dist_cross,
       !is.na(cross_dist_nearest))$cross_dist_nearest))+
     geom_vline(xintercept=0.15,color="grey")
   if(!is.na(threshold)){
       g = g + geom_vline(xintercept=threshold, color="firebrick", linetype=2)
   }
   plots[[sample_id]] = g
   print(g)
}
dev.off()

writeChangeoDb(h,"processed/heavy_chain.tsv")
writeChangeoDb(l,"processed/light_chain.tsv")

# group sequences into clonal clusters
call = paste0("DefineClones.py -d processed/heavy_chain.tsv --dist ",
        0.15," --nproc ",3," --act set --model ham --norm len --gf sample_id")
system(call)

# split clonal clusters based on light chain V/J
call = paste0("~/Programs/immcantation/scripts/light_cluster.py -d processed/heavy_chain_clone-pass.tsv -e processed/light_chain.tsv --doublets count",
            " -o processed/heavy_chain_clone-pass_light-pass.tsv")
system(call)

# check out results
cloned = readChangeoDb("processed/heavy_chain_clone-pass_light-pass.tsv")
max(colSums(table(cloned$sample_id,cloned$clone_id) > 0))
max(table(cloned$clone_id))
writeChangeoDb(cloned,"processed/cloned_heavy_chain.tsv")

# reconstruct germline sequences for each clone
command <- paste("CreateGermlines.py -d", "processed/cloned_heavy_chain.tsv", "-r", r, "--format airr",
	"--cloned --outname cloned_heavy -g dmask")
system(command)

cloned_germ = readChangeoDb("processed/cloned_heavy_germ-pass.tsv")

# calculate SHM frequency in the V gene
cloned <- observedMutations(cloned_germ, sequenceColumn="sequence_alignment",
       germlineColumn="germline_alignment_d_mask",
       regionDefinition=IMGT_V,
       frequency=TRUE,
       combine=TRUE, 
       nproc=3)

writeChangeoDb(cloned,"processed/all_cloned_data.tsv")

