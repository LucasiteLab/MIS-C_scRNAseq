library(Seurat)
library(ggplot2)
library(dplyr)

## Subsetting and scaling 

misc.tcell <- readRDS("shared/t_cell/misc_integrated_tcell_updated_cond.rds")
Idents(misc.tcell) <- "subcluster_annotation1"

misc.nk.cd56dim <- subset(misc.tcell, idents= c("CD56dim_S100A4_NKcells", "CD56dim_CD38NKcells"))
misc.nk.cd56b <- subset(misc.tcell, idents= c("CD56bright_NKcells"))

Idents(misc.nk.cd56dim) <- "condition_new"
Idents(misc.nk.cd56b) <- "condition_new"
misc.nk.cd56dim.peds <- subset(misc.nk.cd56dim, idents = c("C.HD", "MIS-C", "MIS-C-R"))
misc.nk.cd56b.peds <- subset(misc.nk.cd56b, idents = c("C.HD", "MIS-C", "MIS-C-R"))
misc.nk.cd56dim.adult <- subset(misc.nk.cd56dim, idents = c("A.HD", "COVID19-A", "COVID19-B"))
misc.nk.cd56b.adult <- subset(misc.nk.cd56b, idents = c("A.HD", "COVID19-A", "COVID19-B"))

misc.nk.cd56dim.peds <- ScaleData(misc.nk.cd56dim.peds, assay = "RNA")
misc.nk.cd56b.peds <- ScaleData(misc.nk.cd56b.peds, assay = "RNA")
misc.nk.cd56dim.adult <- ScaleData(misc.nk.cd56dim.adult, assay = "RNA")
misc.nk.cd56b.adult <- ScaleData(misc.nk.cd56b.adult, assay = "RNA")

## Make box plot

# CD56dim
Idents(misc.nk.cd56dim.adult) <- "sample_id"
avg.out <- AverageExpression(misc.nk.cd56dim.adult, assays = "RNA", 
                             slot= "scale.data", 
                             features = c("GZMH", "GZMA", "CCL4", "PRF1")) #remember to scale

avg.df <- avg.out$RNA

write.csv(avg.df, file = 
            "Figure_tcell/formatted_plots/GZMH_GZMA_CCL4_PRF1_for_boxplot_CD56dimNK_adult.csv")

# CD56bright
Idents(misc.nk.cd56b.adult) <- "sample_id"
avg.out <- AverageExpression(misc.nk.cd56b.adult, assays = "RNA", 
                             slot= "scale.data", 
                             features = c("GZMH", "GZMA", "CCL4", "PRF1")) #remember to scale

avg.df <- avg.out$RNA

write.csv(avg.df, file = 
            "Figure_tcell/formatted_plots/GZMH_GZMA_CCL4_PRF1_for_boxplot_CD56brightNK_adult.csv")



## Plotting locally peds

library(ggplot2)
library(dplyr)

gene <- read.csv("Sheets/ITGB7_for_boxplot_cd8_memory.csv") #ITGB7 or any gene calculated by AverageExpression 
rownames(gene) <- gene[,1]
gene[,1] <- NULL
gene_t <- as.data.frame(t(as.matrix(gene)))
gene_t$sample <- rownames(gene_t)

# Peds

gene_t[gene_t$sample %in% c("P1.1", "P2.1", "P3.1", "P4.1", "P5.1", "P6.1", "P7.1"), 'condition'] <- 'MIS-C'

gene_t[gene_t$sample %in% c("C.HD1", "C.HD2", "C.HD3", "C.HD4", "C.HD5", "C.HD6"), 'condition'] <- 'C.HD'

gene_t[gene_t$sample %in% c("P3.2", "P4.2"), 'condition'] <- 'MIS-C-R'

gene_severe <- gene_t %>% filter(sample %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1"))


gene_chd <- gene_t %>% filter(condition == "C.HD")
gene_misc <- gene_t %>% filter(condition == "MIS-C")


# Wilcox
w.test <- wilcox.test(x = gene_chd$CCL4, y = gene_misc$CCL4, alternative = c("two.sided"), correct = FALSE)
pval_ccl4 <- w.test$p.value #CCL4 0.03
w.test <- wilcox.test(x = gene_chd$GZMH, y = gene_misc$GZMH, alternative = c("two.sided"), correct = FALSE)
pval_gzmh <- w.test$p.value #GZMH 0.02
w.test <- wilcox.test(x = gene_chd$GZMA, y = gene_misc$GZMA, alternative = c("two.sided"), correct = FALSE)
pval_gzma <- w.test$p.value #GZMA ns
w.test <- wilcox.test(x = gene_chd$PRF1, y = gene_misc$PRF1, alternative = c("two.sided"), correct = FALSE)
pval_prf1 <- w.test$p.value #PRF1 0.005
w.test <- wilcox.test(x = gene_chd$ITGB7, y = gene_misc$ITGB7, alternative = c("two.sided"), correct = FALSE)
pval_itgb7 <- w.test$p.value


# Levels

gene_t$condition <- factor(gene_t$condition, levels = c("C.HD", "MIS-C", "MIS-C-R"))
gene_severe$condition <- factor(gene_severe$condition, levels = c("C.HD", "MIS-C", "MIS-C-R"))
gene_t_2 <- gene_t %>% filter(!(sample %in% c("P1.1", "P2.1", "P3.1", "P6.1", "P7.1")))


cols <- c("#6baed6", "#FC9272", "#969696", "#9970ab", "#ec7014", "#fec44f")

plot1 <- ggplot(gene_t, aes(x = condition, y = ITGB7)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(data=gene_severe, colour ="#c94040",  size = 0.5, width = 0.12)+
  geom_jitter(data = gene_t_2, aes(colour = condition), size = 0.5, width = 0.12) +
  ggtitle("ITGB7") +
  xlab("") +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C")), annotation = "0.02", size = 0.18, textsize = 1.8,
              y_position = 0.8) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line=element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15),
        plot.title = element_text(hjust = 0.5)) +
  ylim(-0.75, 1) +
  ylab("Scaled avg. expression")

plot1

ggsave(plot1, file = "New_formatting/ITGB7_CD8mem_peds_new.pdf", height = 1.5, width =1)
