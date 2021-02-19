library(reshape2)
library(tidyverse)
library(ggplot2)

samples <- read.csv("Sheets/sample_myeloid_compositions.txt", sep = "\t", row.names = NULL)

for(j in 2:ncol(samples)){
  samples[10,j] <- print(sum(samples[1:9,j]))
}

samples_prop <- samples

#divide by donor contributuion
for(j in 2:ncol(samples)){ 
  for(i in 1:nrow(samples)){
    samples_prop[i,j] <- print(samples[i,j]/samples[nrow(samples),j])
  }
}

names(samples_prop)[1] <- "annotation_myeloid"
donor.names <- colnames(samples_prop)[2:ncol(samples_prop)]
samples_prop <- samples_prop[1:9,]
samples_plot <- melt(samples_prop[1:9,], measure.vars = donor.names, id.vars = c("annotation_myeloid"))

# Assign conditions

samples_plot[samples_plot$variable %in% c("C27", "C32", "C33", "C39", "HA5876", "HA5877",
                                          "HA5894", "HA5952", "HA5953", "HA5957", "HD_32M", 
                                          "HD_35F", "HD_36M"), 'condition'] <- 'A.HD'

samples_plot[samples_plot$variable %in% c("NS0A", "NS1A", "TS2A","TS3A"), 'condition'] <- 'COVID19-A'

samples_plot[samples_plot$variable %in% c("NS0B", "NS1B", "TP8B", "TP9B", 
                                          "TS2B", "TS3B"), 'condition'] <- 'COVID19-B'

samples_plot[samples_plot$variable %in% c("Y124.1", "Y125.1", "Y127.1",
                                          "Y111.1", "Y113.1", "Y117.1", "Y129.1"), 'condition'] <- 'MIS-C'

samples_plot[samples_plot$variable %in% c("NC.13F", "Y28.2", "Y28.4", 
                                          "Y29.2", "Y54.4", "Y70.4"), 'condition'] <- 'C.HD'

samples_plot[samples_plot$variable %in% c("Y117.R", "Y124.R"), 'condition'] <- 'MIS-C-R'


samples_plot_severe <- samples_plot %>% filter(variable %in% c("Y111.1", "Y113.1", "Y117.1", "Y127.1","Y129.1"))

samples_plot[,5] <- samples_plot[,3]*100 #change to 6 when we have annotations
samples_plot_severe[,5] <- samples_plot_severe[,3]*100 #change to 6 when we have annotations
names(samples_plot)[5] <- "percentage"
names(samples_plot_severe)[5] <- "percentage"

level_order <-  c("C.HD", "MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")

annotation_order <- samples_prop[,1]

samples_plot$annotation_myeloid <- factor(samples_plot$annotation_myeloid, level = annotation_order)
samples_plot_severe$annotation_myeloid <- factor(samples_plot_severe$annotation_myeloid, level = annotation_order)

samples_plot <- samples_plot %>% filter(condition %in% c("C.HD", "MIS-C", "MIS-C-R"))
samples_plot_severe <- samples_plot_severe %>% filter(condition %in% c("C.HD", "MIS-C", "MIS-C-R"))

##

condition_pvals <- list()

pvals <- list()

level_order <- c("C.HD", "MIS-C", "MIS-C-R")

cluster_num <- annotation_order

for(j in 1:length(level_order)){
  for(i in 1:length(cluster_num)){
    cluster <- samples_plot %>% filter(annotation_myeloid == cluster_num[i] & condition == level_order[j])
    pvals[[i]] <- cluster$percentage
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

names(chd_misc_df)[2] <- "annotation_myeloid"

#write.csv(chd_misc_df2, file = "wilcoxon_pbmc_compositions.csv")

# 
chd_misc_df$comp.a <- 'C.HD'
chd_misc_df$comp.b <- 'MIS-C'
#samples_plot2 <- samples_plot
#samples_plot2$condition <- as.character(samples_plot2$condition)
chd_misc_df$ypos <- c(35,60,32,25,25,20,22,12,10)
chd_misc_df$annotation_myeloid <- as.character(chd_misc_df$annotation_myeloid)
#chd_misc_df$pval <- format.pval(as.double(as.character(chd_misc_df2$pval)), digitsc=1, eps=0.001)
chd_misc_df$annotation_myeloid <- factor(chd_misc_df$annotation_myeloid, level = annotation_order)



for(i in 1:nrow(chd_misc_df)){
  pval <- chd_misc_df[i,1]
  chd_misc_df[i,1] <- format.pval(as.double(as.character(pval)), digits=1, eps=0.001)
}


chd_misc_df <- chd_misc_df %>% filter(pval < 0.051)
samples_plot2 <- samples_plot %>% filter(!(variable %in% c("Y111.1", "Y113.1", "Y117.1", "Y127.1","Y129.1")))

cols <- c("#6baed6", "#FC9272", "#969696")

plot1 <- ggplot(samples_plot, aes(x = factor(condition, level = level_order), y = percentage)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(data=samples_plot2, aes(colour = factor(condition, level = level_order)), size = 3, width = 0.25) + 
  geom_jitter(data=samples_plot_severe, aes(group = factor(condition, level = level_order)), colour ="#c94040",  size = 3, width = 0.25)+
  geom_signif(data=chd_misc_df, aes(xmin=comp.a, xmax=comp.b, annotations=pval, y_position=ypos),
              textsize = 5, vjust=-0.2, manual=TRUE) + ylim(c(NA, 5)) +
  facet_wrap(~annotation_myeloid, nrow=2, scales='free_y', labeller = label_wrap_gen(16)) +
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
        strip.text = element_text(size = 13, margin =margin()),
        strip.text.x = element_text(margin = margin( b = 10, t = 10)),
        strip.background = element_rect(fill="white"),
        legend.position = c(0.9,0.07)) +
  ylab("Cell Percentage") +
  xlab("") 
plot1

ggsave(plot = plot1, "Sample_composition_graphs_resubmit/sample_composition_grid_myeloid_severe_mod_final.pdf", height= 7, width = 10)
