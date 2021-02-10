library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsignf)

misc.myeloid <- readRDS("shared/Myeloid/misc_integrated_myeloid_updated_cond.rds")

# Subsetting on population of interest
misc.neut <- subset(misc.myeloid, idents = c("Neutrophils I", "Neutrophils II"))
misc.mono <- subset(misc.myeloid, idents= c("Classical Monocytes I", "Classical Monocytes II", 
                                            "Classical Monocytes III", "Intermediate Monocytes",
                                            "Non Classical Monocytes"))
misc.c.mono <- subset(misc.myeloid, idents= c("Classical Monocytes I", "Classical Monocytes II", 
                                              "Classical Monocytes III"))
misc.nc.int.mono <- subset(misc.myeloid, idents= c("Intermediate Monocytes","Non Classical Monocytes"))

misc.sub.myeloid <- subset(misc.myeloid, idents = c("Classical Monocytes I", "Classical Monocytes II", 
                                                    "Classical Monocytes III", "Intermediate Monocytes",
                                                    "Non Classical Monocytes", "Neutrophils I",
                                                    "Neutrophils II"))


# Alarmin Module score

alarmin <- list(c("S100A8", "S100A9", "S100A12"))

misc.cluster <- AddModuleScore(misc.sub.myeloid, name = "Alarmin_score", nbins=24, ctrl=100,
                               features = alarmin, assay = "RNA") #not scaled

cols <- c("#6baed6", "#c94040", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

condition_order <-  c("C.HD", "MIS-C-S", "MIS-C-M", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")
misc.cluster@meta.data$condition_severe <- factor(misc.cluster@meta.data$condition_severe, level = condition_order)
sample_order <-  c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6",
                   "P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1",
                   "P3.2", "P4.2", "A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                   "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                   "A.HD13", "A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1", "A.COV1.2",
                   "A.COV2.2", "A.COV3.2", "A.COV4.2", "A.COV5.2", "A.COV6.2")
misc.cluster@meta.data$sample_id <- factor(misc.cluster@meta.data$sample_id, level = sample_order)
#pbmc_meta <- misc.cluster[[]]

Idents(misc.cluster) <- "condition_new"
misc.peds <- subset(misc.cluster, idents = c("C.HD", "MIS-C", "MIS-C-R"))
misc.adult <- subset(misc.cluster, idents = c("A.HD", "COVID19-A", "COVID19-B"))

# Violin plot (not pictured in figure)

bacplot <- VlnPlot(misc.peds, group.by = "sample_id", split.by = "condition_severe",
                   features = "Alarmin_score1", pt.size = 0) +
  ggtitle("S100A8.A9.A12") + 
  scale_fill_manual(values = cols) + theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ylab("Module score") +
  xlab(" ") 

q <- ggplot_build(bacplot)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.2
#q$data[[2]]$alpha <- 0.5
bacplot <- ggplot_gtable(q)


ggsave(plot = bacplot, "Figure_myeloid/formatted_figures/alarmin_cmono_cohort_peds_severe_new.pdf", height = 1.5, width = 2)

meta <- misc.adult[[]]
df_score <- data.frame('sample' = meta$sample_id, 'score' = meta$Alarmin_score1, 
                       'cluster' = meta$annotation, 'condition' = meta$condition_new)
write.csv(df_score, file = "Figure_myeloid/alarmin_myeloid_specific_module_adult.csv") # Export to plot


# HLA score

hla <- list(c("HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2",
              "HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5"))

misc.cluster <- AddModuleScore(misc.mono, name = "HLA_score", nbins=24, ctrl=100,
                               features = hla, assay = "RNA") #not scaled

cols <- c("#6baed6", "#FC9272", "#c94040", "#969696", "#9970ab", "#ec7014", "#fec44f")
cols2 <- c("#9970ab", "#ec7014", "#fec44f")

condition_order <-  c("C.HD", "MIS-C-S", "MIS-C-M", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

misc.cluster@meta.data$condition_severe <- factor(misc.cluster@meta.data$condition_severe, level = condition_order)

sample_order <-  c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6",
                   "P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1",
                   "P3.2", "P4.2", "A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                   "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                   "A.HD13", "A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1", "A.COV1.2",
                   "A.COV2.2", "A.COV3.2", "A.COV4.2", "A.COV5.2", "A.COV6.2")

misc.cluster@meta.data$sample_id <- factor(misc.cluster@meta.data$sample_id, level = sample_order)
#pbmc_meta <- misc.cluster[[]]

Idents(misc.cluster) <- "condition_new"
misc.peds <- subset(misc.cluster, idents = c("C.HD", "MIS-C", "MIS-C-R"))
misc.adult <- subset(misc.cluster, idents = c("A.HD", "COVID19-A", "COVID19-B"))

# If by seurat violin plot
bacplot <- VlnPlot(misc.peds, group.by = "sample_id", split.by = "condition_new",
                   features = "HLA_score1", pt.size = 0) +
  ggtitle("HLA class II score in monocytes") + 
  scale_fill_manual(values = cols) + theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ylab("Scaled expression level") +
  xlab(" ") 

q <- ggplot_build(bacplot)
#q$data[[1]]$colour <- NA
q$data[[1]]$size <- 0.2
#q$data[[2]]$alpha <- 0.5
bacplot <- ggplot_gtable(q)

ggsave(plot = bacplot, "Figure_myeloid/formatted_figures/hla_neut_peds_donor.pdf", height = 1.5, width = 2)


meta <- misc.peds[[]]
df_score <- data.frame('sample' = meta$sample_id, 'score' = meta$HLA_score1, 
                       'cluster' = meta$annotation, 'condition' = meta$condition_new)
write.csv(df_score, file = "Figure_myeloid/hla_myeloid_module_peds_neut_mono.csv")

saveRDS(misc.peds, "monocyte_pediatric_alarmin_score.rds")

# Individual genes

# Monocytes
Idents(misc.mono) <- "condition_new"
misc.peds <- subset(misc.mono, idents = c("C.HD", "MIS-C", "MIS-C-R"))
misc.adult <- subset(misc.mono, idents = c("A.HD", "COVID19-A", "COVID19-B"))
misc.peds <- ScaleData(misc.peds, assay = "RNA")
misc.adult <- ScaleData(misc.adult, assay = "RNA")

condition_order <-  c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

misc.adult@meta.data$condition_new <- factor(misc.adult@meta.data$condition_new, level = condition_order)

sample_order <-  c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6",
                   "P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1",
                   "P3.2", "P4.2", "A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                   "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                   "A.HD13", "A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1", "A.COV1.2",
                   "A.COV2.2", "A.COV3.2", "A.COV4.2", "A.COV5.2", "A.COV6.2")

misc.adult@meta.data$sample_id <- factor(misc.adult@meta.data$sample_id, level = sample_order)

saveRDS(misc.adult, "Figure_myeloid/formatted_figures/myeloid_adult_obj_scaled.rds")


# Making table
Idents(misc.adult) <- "sample_id"
avg.out <- AverageExpression(misc.adult, assays = "RNA", 
                             slot= "scale.data", 
                             features = c("CD86")) #remember to scale

avg.df <- avg.out$RNA

write.csv(avg.df, file = "Figure_myeloid/cd86_monocyte_adult.csv")


## Plotting individual genes locally

gene <- read.csv("Sheets/cd86_monocyte_adult.csv") #gene just a placeholder name
rownames(gene) <- gene[,1]
gene[,1] <- NULL
gene_t <- as.data.frame(t(as.matrix(gene)))
gene_t$sample <- rownames(gene_t)

# Adults
gene_t[gene_t$sample %in% c("A.HD1", "A.HD2", "A.HD3", "A.HD4", "A.HD5",
                            "A.HD6", "A.HD7", "A.HD8", "A.HD9", "A.HD10", "A.HD11", "A.HD12", 
                            "A.HD13"), 'condition'] <- 'A.HD'

gene_t[gene_t$sample %in% c("A.COV1.1", "A.COV2.1", "A.COV3.1", "A.COV4.1"), 'condition'] <- 'COVID19-A'

gene_t[gene_t$sample %in% c("A.COV1.2","A.COV2.2", "A.COV3.2", 
                            "A.COV4.2", "A.COV5.2", "A.COV6.2"), 'condition'] <- 'COVID19-B'

gene_ahd <- gene_t %>% filter(condition == "A.HD")
gene_covida <- gene_t %>% filter(condition == "COVID19-A")
gene_covidb <- gene_t %>% filter(condition == "COVID19-B")

# Peds

gene_t[gene_t$sample %in% c("P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1"), 'condition'] <- 'MIS-C'

gene_t[gene_t$sample %in% c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6"), 'condition'] <- 'C.HD'

gene_t[gene_t$sample %in% c("P3.2", "P4.2"), 'condition'] <- 'MIS-C-R'

gene_severe <- gene_t %>% filter(sample %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1"))

gene_chd <- gene_t %>% filter(condition == "C.HD")
gene_misc <- gene_t %>% filter(condition == "MIS-C")


# Wilcox
w.test <- wilcox.test(x = gene_ahd$CD86, y = gene_covidb$CD86, alternative = c("two.sided"), correct = FALSE)
pval_covidb <- w.test$p.value

w.test <- wilcox.test(x = gene_chd$CD86, y = gene_misc$CD86, alternative = c("two.sided"), correct = FALSE)
pval_cd86 <- w.test$p.value


# Levels

gene_t$condition <- factor(gene_t$condition, levels = c("A.HD", "COVID19-A", "COVID19-B"))


gene_t$condition <- factor(gene_t$condition, levels = c("C.HD", "MIS-C", "MIS-C-R"))
gene_severe$condition <- factor(gene_severe$condition, levels = c("C.HD", "MIS-C", "MIS-C-R"))


cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")
cols2 <- c("#9970ab", "#ec7014", "#fec44f")

gene_t_2 <- gene_t %>% filter(!(sample %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))


plot1 <- ggplot(gene_t, aes(x = condition, y = CD86)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(data=gene_severe, colour ="#c94040",  size = 0.5, width = 0.12) +
  geom_jitter(data = gene_t_2, aes(colour = condition), size = 0.5, width = 0.12) +
  ggtitle("CD86") +
  xlab("") +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C")), annotation = "0.002", 
              size = 0.12, textsize = 2,
              y_position = 1)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylim(-0.75, 1.5) +
  ylab("Scaled avg. expression")
plot1

ggsave(plot1, file = "New_formatting/CD86_peds_monocyte.pdf", height = 1.5, width =1)


#Adult
plot1 <- ggplot(gene_t, aes(x = condition, y = CD86)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(aes(colour = condition), size = 0.5, width = 0.12) +
  ggtitle("CD86") +
  xlab("") +
  scale_color_manual(values = cols2) +
  geom_signif(comparisons = list(c("A.HD", "COVID19-A")), annotation = "ns", size = 0.12, 
              textsize = 2, y_position = 0.5)+
  geom_signif(comparisons = list(c("A.HD", "COVID19-B")), annotation = "0.02", size = 0.12, 
              textsize = 2, y_position = 0.75)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylim(-0.5, 1) +
  ylab("Scaled avg. expression")
plot1

ggsave(plot1, file = "New_formatting/CD86_adult_monocyte.pdf", height = 1.5, width =1)



## Signature plotting peds

##
cyto <- read.csv("Sheets/hla_monocyte_module_peds.csv")

names <- unique(cyto$sample)
means <- data.frame("name" = rep(NA,38), "value" = rep(NA,38))

for(i in 1:length(names)){
  cyto_tmp <- cyto %>% filter(sample == names[i])
  means[i,1] <- names[i]
  means[i,2] <- mean(cyto_tmp[,3])
}


means[means$name %in% c("P1.1", "P2.1", "P3.1", "P4.1",
                        "P5.1", "P6.1", "P7.1"), 'condition'] <- 'MIS-C'

means[means$name %in% c("C.HD1", "C.HD2", "C.HD3", 
                        "C.HD4", "C.HD5", "C.HD6"), 'condition'] <- 'C.HD'

means[means$name %in% c("P3.2", "P4.2"), 'condition'] <- 'MIS-C-R'

means_severe <- means %>% filter(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1"))


means <- na.omit(means)

level_order <-  c("C.HD", "MIS-C", "MIS-C-R")

means$condition <- factor(means$condition, level = level_order)
means_severe$condition <- factor(means_severe$condition, level = level_order)

misc_tmp <- means %>% filter(condition == "MIS-C")
chd_tmp <- means %>% filter(condition == "C.HD")


# Example

w.test <- wilcox.test(x = chd_tmp$value, y = misc_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval <- w.test$p.value


cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

means2 <- means %>% filter(!(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))

plot1 <- ggplot(means, aes(x = condition, y = value)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(data=means_severe, colour ="#c94040",  size = 0.5, width = 0.12)+
  geom_jitter(data = means2, aes(colour = factor(condition, level = level_order)), 
              size = 0.5, width = 0.12) +
  ggtitle("Class II HLA") +
  xlab("") +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C")), annotation = "0.008", 
              size = 0.12, textsize = 2,
              y_position =1.5)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15), ) +
  ylim(c(-1, 2)) +
  ylab("Average module score")

plot1

ggsave(plot1, file = "New_formatting/mono_hla_score_peds_severe_new.pdf", 
       height = 1.5, width =1)



## Signature plotting adult

##
cyto <- read.csv("Sheets/hla_monocyte_module_adult.csv")

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


means <- na.omit(means)

level_order <-  c("A.HD", "COVID19-A", "COVID19-B")

means$condition <- factor(means$condition, level = level_order)

covid19a_tmp <- means %>% filter(condition == "COVID19-A")
covid19b_tmp <- means %>% filter(condition == "COVID19-B")
ahd_tmp <- means %>% filter(condition == "A.HD")

# Example

w.test <- wilcox.test(x = ahd_tmp$value, y = covid19a_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval_covida <- w.test$p.value

cols2 <- c("#9970ab", "#ec7014", "#fec44f")


plot1 <- ggplot(means, aes(x = condition, y = value)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(aes(colour = factor(condition, level = level_order)), size = 0.5, width = 0.12) +
  ggtitle("Class II HLA") +
  xlab("") +
  scale_color_manual(values = cols2) +
  geom_signif(comparisons = list(c("A.HD", "COVID19-A")), annotation = "ns", size = 0.12, textsize = 2,
              y_position = 1.3)+
  geom_signif(comparisons = list(c("A.HD", "COVID19-B")), annotation = "0.007", size = 0.12, textsize = 2,
              y_position = 1.7)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylim(c(-0.75, 2)) +
  ylab("Average module score")

plot1

ggsave(plot1, file = "New_formatting/mono_hla_score_adult_jittered.pdf", 
       height = 1.5, width =1)
