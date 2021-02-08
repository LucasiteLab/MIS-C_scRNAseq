library(Seurat)
library(tidyverse)
library(ggplot2)
library(AUCell)

## Module score for Viral across PBMCs and myeloid cells (Tirosh et al. Science 2016)

misc.cluster <- readRDS("shared/pbmc_2.0/misc_integrated_object_2.0_updated_cond.rds")

# Modify factor levels
condition_order <-  c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")
misc.cluster@meta.data$condition_new <- factor(misc.cluster@meta.data$condition_new, level = condition_order)
sample_order <-  c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6",
                   "P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1",
                   "P3.2", "P4.2", "A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                   "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                   "A.HD13", "A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1", "A.COV1.2",
                   "A.COV2.2", "A.COV3.2", "A.COV4.2", "A.COV5.2", "A.COV6.2")
misc.cluster@meta.data$sample_id <- factor(misc.cluster@meta.data$sample_id, level = sample_order)
cols <- c("#6baed6", "#c94040", "#969696", "#9970ab", "#ec7014", "#fec44f")

# Viral and bacterial scores from Lydon et al. (Respiratory infections)
ifn <- list(c("SIGLEC1","RABGAP1L","IFI27","CADM1","RSAD2","MX1","SERPING1")) 
bac <- list(c("SMPD1","CD44","SERPING1","SPI1","HERC1","MCTP1","FOLR3","CFAP45",
              "PRF1","CTBP1","HLA-DRB1","ARL1","OAS3","ZER1", "IFIT2","IFITM1"))
viral_h <- list(c("IFI44L", "IFI27", "RSAD2", "SIGLEC1", "IFIT1", "ISG15"))
bac2 <- list(c("HK3", "TNIP1", "GPAA1", "CTSB"))
sepsis <- list(c("PLAC8", "CLU", "RETN", "CD63", "ALOX5AP", "SEC61G", "TXN", "MT1X"))

# Identify monocyte clusters

misc.cluster <- AddModuleScore(misc.cluster, name = "Viral_score_up", nbins=24, ctrl=100,
                             features = sepsis, assay = "RNA") #not scaled

Idents(misc.cluster) <- "seurat_clusters"
misc.mono <- subset(misc.cluster, idents = c(2,8,12)) #just mono and neut
misc.mono <- subset(misc.cluster, idents = c(2,12)) #just mono

pbmc_meta <- misc.mono[[]]

# Export for box plot statistical analysis

df_myeloid <- data.frame('sample' = pbmc_meta$sample_id, 'score' = pbmc_meta$Viral_score_up1, 
                         'cluster' = pbmc_meta$annotation, 'condition' = pbmc_meta$condition)
write.csv(df_myeloid, file = "AUCell/sepsis_mono_module_score.csv") #use to calculate pvalue 


#
library(dplyr)
library(ggplot2)

cyto <- read.csv("Sheets/sepsis_mono_module_score.csv")

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

w.test <- wilcox.test(x = covid19a_tmp$value, y = misc_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval <- w.test$p.value


means2 <- means %>% filter(!(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))

cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

plot1 <- ggplot(means, aes(x = condition, y = value)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) + 
  geom_jitter(data=means_severe, colour ="#c94040",  size = 0.5, width = 0.1)+
  geom_jitter(data = means2, aes(colour = factor(condition, level = level_order)), size = 0.5, width = 0.1) +
  scale_color_manual(values = cols) +
  ggtitle("Sepsis associated monocyte sig") +
  xlab("") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 8, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylab("Average module score") 

plot1

ggsave(plot1, file = "Fig2_viral_bac_fig/Sepsis_monocyte_module_Hacohen.pdf", height = 2, width =3)







## AUCell


