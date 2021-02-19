library(reshape2)
library(Hmisc)
library(stats)
library(gtools)
library(tidyverse)


#  Cell frequencies


samples <- read.csv("Sheets/sample_full_integrated_pbmc_wo112.csv")

for(j in 2:ncol(samples)){
  samples[31,j] <- print(sum(samples[1:30,j]))
}

samples_prop <- samples


for(j in 2:ncol(samples)){ 
  for(i in 1:nrow(samples)){
    samples_prop[i,j] <- print(samples[i,j]/samples[nrow(samples),j])
  }
}

names(samples_prop)[1] <- "seurat_clusters"
donor.names <- colnames(samples_prop)[2:ncol(samples_prop)]
samples_prop <- samples_prop[1:30,]


annotations <- read.csv("Sheets/working_annotation_2.0.csv")
samples_new <- cbind(samples_prop[,2:ncol(samples_prop)], annotations=annotations[,2])
rownames(samples_new) <- samples_new$annotations
samples_new[,39] <- NULL

clin_s <- t(as.matrix(samples_new))
clin_m <- subset(clin_s, rownames(clin_s) %in% c("Y111.1", "Y113.1", "Y117.1", "Y124.1",
                                                 "Y125.1", "Y127.1", "Y129.1"))

clin_c <- subset(clin_s, rownames(clin_s) %in% c("NC.13F", "Y28.2", "Y28.4", 
                                                 "Y29.2", "Y54.4", "Y70.4"))

torem <- c(19,22,27,28)
#clin_s <- as.data.frame(clin_s)
clin_c <- as.data.frame(clin_c)
clin_m <- as.data.frame(clin_m)

#clin_s[,torem] <- NULL

clin_c[,torem] <- NULL
clin_m[,torem] <- NULL
#clin_s <- as.matrix(clin_s)
clin_c <- as.matrix(clin_c)
clin_m <- as.matrix(clin_m)

#cormatrix <- rcorr(clin_s, type = "spearman")
cormatrix_c <- rcorr(clin_c, type = "spearman")
cormatrix_m <- rcorr(clin_m, type = "spearman")

#padj <- matrix(p.adjust(cormatrix$P, method='BH'), nrow=nrow(cormatrix$P), ncol=ncol(cormatrix$P)) # for existing
padj_c <- matrix(p.adjust(cormatrix_c$P, method='BH'), nrow=nrow(cormatrix_c$P), ncol=ncol(cormatrix_c$P)) 
padj_m <- matrix(p.adjust(cormatrix_m$P, method='BH'), nrow=nrow(cormatrix_m$P), ncol=ncol(cormatrix_m$P)) 

cormatrix_c$P.adj <- padj_c
cormatrix_m$P.adj <- padj_m

col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

pdf(file = "corrplot_cell_freq_chd_only.pdf") #if using clin_m, then MIS-C samples only
corrplot(cormatrix_c$r, method = "square",  order = "hclust", tl.col = "black", 
         col = col2(50), p.mat = cormatrix_c$P.adj, insig = "label_sig", 
         sig.level = c(.001, .01, .05),
         pch.cex = .5, pch.col = "white", diag = FALSE)
dev.off()

