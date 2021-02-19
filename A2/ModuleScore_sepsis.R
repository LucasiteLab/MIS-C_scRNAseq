library(Seurat)
library(tidyverse)
library(ggplot2)

## Module score for sepsis across classical monocytes

misc.cluster <- readRDS("shared/Myeloid/misc_integrated_myeloid_updated_cond.rds")

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
sepsis <- list(c("PLAC8", "CLU", "RETN", "CD63", "ALOX5AP", "SEC61G", "TXN", "MT1X"))


# Identify monocyte clusters

misc.cluster <- AddModuleScore(misc.cluster, name = "Viral_score_up", nbins=24, ctrl=100,
                               features = sepsis, assay = "RNA") #not scaled


misc.mono <- subset(misc.cluster, idents = c("Classical Monocytes I",
                                             "Classical Monocytes II", "Classical Monocytes III")) #just mono and neut


pbmc_meta <- misc.mono[[]]

# Export for box plot statistical analysis

df_myeloid <- data.frame('sample' = pbmc_meta$sample_id, 'score' = pbmc_meta$Viral_score_up1, 
                         'cluster' = pbmc_meta$annotation, 'condition' = pbmc_meta$condition)
write.csv(df_myeloid, file = "ModuleScore/sepsis_module_score_cmono_subclusters.csv") #use to calculate pvalue 


#
library(dplyr)
library(ggplot2)
library(ggsignif)

cyto <- read.csv("Sheets/sepsis_module_score_cmono_subclusters.csv")

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
pval_ped <- w.test$p.value
w.test <- wilcox.test(x = ahd_tmp$value, y = covid19a_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval_covida <- w.test$p.value
w.test <- wilcox.test(x = ahd_tmp$value, y = covid19b_tmp$value, alternative = c("two.sided"), correct = FALSE)
pval_covidb <- w.test$p.value


means2 <- means %>% filter(!(name %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))

cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

plot1 <- ggplot(means, aes(x = condition, y = value)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) + 
  geom_jitter(data=means_severe, colour ="#c94040",  size = 0.5, width = 0.1)+
  geom_jitter(data = means2, aes(colour = factor(condition, level = level_order)), size = 0.5, width = 0.1) +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C")), annotation = "0.03", 
              size = 0.12, textsize = 2,
              y_position =0.2)+
  geom_signif(comparisons = list(c("A.HD", "COVID19-A")), annotation = "0.0008", 
              size = 0.12, textsize = 2,
              y_position =0.32)+
  geom_signif(comparisons = list(c("A.HD", "COVID19-B")), annotation = "<0.0001", 
              size = 0.12, textsize = 2,
              y_position =0.45)+
  ggtitle("Sepsis signature in classical monocytes") +
  xlab("") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 7, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylab("Average module score") +
  ylim(c(-0.6, 0.5))

plot1

ggsave(plot1, file = "sepsis_monocyte_module_score_subclusters.pdf", height = 2, width =3)
