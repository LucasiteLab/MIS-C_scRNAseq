library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggsignif)

samples <- read.csv("Sheets/cell_numbers_NKT_formatted.csv")

rownames(samples) <- samples[,1]
samples[,1] <- NULL

samples_t  <- as.data.frame(t(as.matrix(samples)))
samples_t$annotation_tcell <- rownames(samples_t)

samples <- samples_t[, c(39, 1:38)]

samples <- samples %>% filter(annotation_tcell != "T.NK.monocyte_doublets")

samples$annotation_tcell <- gsub("_", " ", samples$annotation_tcell)

for(j in 2:ncol(samples)){
  samples[18,j] <- print(sum(samples[1:17,j]))
}

samples_prop <- samples

#divide by donor contributuion
for(j in 2:ncol(samples)){ 
  for(i in 1:nrow(samples)){
    samples_prop[i,j] <- print(samples[i,j]/samples[nrow(samples),j])
  }
}

#names(samples_prop)[1] <- "annotation_myeloid"
donor.names <- colnames(samples_prop)[2:ncol(samples_prop)]
samples_prop <- samples_prop[1:17,]
samples_plot <- melt(samples_prop[1:17,], measure.vars = donor.names)



samples_plot[samples_plot$variable %in% c("C27", "C32", "C33", "C39", "HA5876", "HA5877",
                                              "HA5894", "HA5952", "HA5953", "HA5957", "HD_32M", 
                                              "HD_35F", "HD_36M"), 'condition'] <- 'A.HD'

samples_plot[samples_plot$variable %in% c("NS0A", "NS1A", "TS2A","TS3A"), 'condition'] <- 'COVID19-A'

samples_plot[samples_plot$variable %in% c("NS0B", "NS1B", "TP8B", "TP9B", 
                                              "TS2B", "TS3B"), 'condition'] <- 'COVID19-B'

samples_plot[samples_plot$variable %in% c("Y124-1", "Y125-1", "Y127-1",
                                              "Y111-1", "Y113-1", "Y117-1", "Y129-1"), 'condition'] <- 'MIS-C'

samples_plot[samples_plot$variable %in% c("NC-13F", "Y28-2", "Y28-4", 
                                              "Y29-2", "Y54-4", "Y70-4"), 'condition'] <- 'C.HD'

samples_plot[samples_plot$variable %in% c("Y117-R", "Y124-R"), 'condition'] <- 'MIS-C-R'

samples_plot_severe <- samples_plot %>% filter(variable %in% c("Y111-1", "Y113-1", "Y117-1", "Y127-1","Y129-1"))

samples_plot[,5] <- samples_plot[,3]*100 #change to 6 when we have annotations
samples_plot_severe[,5] <- samples_plot_severe[,3]*100 #change to 6 when we have annotations

names(samples_plot)[5] <- "percentage"
names(samples_plot_severe)[5] <- "percentage"


#


level_order <-  c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

annotation_order <- c("Naive CD4 Tcells 1", "Naive CD4 Tcells 2", "Naive CD8 Tcells", "Naive CD4 Tcells 3", 
                      "Memory CD4 Tcells", "CXCR3 Memory CD4 Tcells", "CCR6 Memory CD4 Tcells", 
                      "Effector memory CD8 Tcells", "Central memory CD8 Tcells", "Terminal effector memory CD8 Tcells", 
                      "regulatory Tcells", "CD56dim CD38NKcells", "CD56bright NKcells", "CD56dim S100A4 NKcells", 
                      "Ki67 NK.Tcells", "Vd2 gdTcells", "MAIT cells")

samples_plot$annotation_tcell <- factor(samples_plot$annotation_tcell, level = annotation_order)
samples_plot_severe$annotation_tcell <- factor(samples_plot_severe$annotation_tcell, level = annotation_order)


colnames(samples_plot) <- c("annotation_tcell", "donor", "value", "condition", "percentage")
colnames(samples_plot_severe) <- c("annotation_tcell", "donor", "value", "condition", "percentage")


samples_plot <- samples_plot %>% filter(condition %in% c("C.HD", "MIS-C", "MIS-C-R"))
samples_plot_severe <- samples_plot_severe %>% filter(condition %in% c("C.HD", "MIS-C", "MIS-C-R"))

###


condition_pvals <- list()

pvals <- list()

cluster_num <- annotation_order

level_order <- c("C.HD", "MIS-C", "MIS-C-R")

for(j in 1:length(level_order)){
  for(i in 1:length(cluster_num)){
    cluster <- samples_plot %>% filter(annotation_tcell == cluster_num[i] & condition == level_order[j])
    pvals[[i]] <- cluster$value
  }
  names(pvals) <- annotation_order
  condition_pvals[[j]] <- pvals
  names(condition_pvals)[j] <- level_order[j]
}


# For CHD/MISC comparison 

chd_misc_pvals <- c()
annot <- names(pvals)

for(i in 1:length(annot)){
  w.test <- wilcox.test(x = condition_pvals[['C.HD']][[annot[i]]],
                        y = condition_pvals[['MIS-C']][[annot[i]]],
                        alternative = c("two.sided"), correct = FALSE)
  chd_misc_pvals[i] <- w.test$p.value
}

chd_misc_df <- data.frame(pval=chd_misc_pvals, facet_name=names(pvals))

names(chd_misc_df)[2] <- "annotation_tcell"

#write.csv(chd_misc_df2, file = "wilcoxon_pbmc_compositions.csv")

# 
chd_misc_df$comp.a <- 'C.HD'
chd_misc_df$comp.b <- 'MIS-C'
#samples_plot2 <- samples_plot
#samples_plot2$condition <- as.character(samples_plot2$condition)
chd_misc_df$ypos <- c(30,18,20,10,8,3,6,6,5, 10, 10,12,4,5,12,8,7)
chd_misc_df$annotation_tcell <- as.character(chd_misc_df$annotation_tcell)
#chd_misc_df$pval <- format.pval(as.double(as.character(chd_misc_df2$pval)), digitsc=1, eps=0.001)
chd_misc_df$annotation_tcell <- factor(chd_misc_df$annotation_tcell, level = annotation_order)



for(i in 1:nrow(chd_misc_df)){
  pval <- chd_misc_df[i,1]
  chd_misc_df[i,1] <- format.pval(as.double(as.character(pval)), digits=1, eps=0.001)
}


chd_misc_df <- chd_misc_df %>% filter(pval < 0.051)

#
samples_plot2 <- samples_plot %>% filter(!(donor %in% c("Y111-1", "Y113-1", "Y117-1", "Y127-1","Y129-1")))
cols <- c("#6baed6", "#FC9272", "#969696")

plot1 <- ggplot(samples_plot, aes(x = factor(condition, level = level_order), y = percentage)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(data= samples_plot2, aes(colour = factor(condition, level = level_order)), size = 3, width = 0.25) + 
  geom_jitter(data=samples_plot_severe, aes(group = factor(condition, level = level_order)), colour ="#c94040",  size = 3, width = 0.25)+
  geom_signif(data=chd_misc_df, aes(xmin=comp.a, xmax=comp.b, annotations=pval, y_position=ypos),
              textsize = 5, vjust=-0.2, manual=TRUE) + ylim(c(NA, 5)) +
  facet_wrap(~annotation_tcell, nrow=3, scales='free_y', labeller = label_wrap_gen(16)) +
  scale_color_manual(values = cols) +theme_bw() +
  scale_y_continuous(expand = expansion(mult=c(0.1, 0.25))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size =10),
        legend.text = element_text(size=16),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, margin =margin()),
        strip.text.x = element_text(margin = margin( b = 3, t = 3)),
        strip.background = element_rect(fill="white"),
        legend.position = c(0.9,0.06)) +
  ylab("Cell Percentage") +
  xlab("") 
plot1

ggsave(plot = plot1, "Sample_composition_graphs_resubmit/sample_composition_grid_tcell_severe_mod_final.pdf", height= 10, width = 14)
