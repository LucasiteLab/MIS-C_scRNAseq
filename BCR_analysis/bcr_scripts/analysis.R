# Kenneth B. Hoehn
# 2/1/2021
# Make all the BCR figures!

library(alakazam)
library(dplyr)
library(tidyr)
library(ggpubr)

cols_sm = c("MIS-C-S"="#C94040", "MIS-C-M"="#FC9272", 
  "C.HD"="#6BAED6", "MIS-C-R"="#969696")

cols <- c("MIS-C"="#C94040","MIS-C-R"="#969696", "C.HD"="#6BAED6",
	"A.HD"="#9970AB", "COVID-A"="#EC7014", "COVID-B"="#FEC44F",
	"ICU"="#969696","MIS-C-S"="#C94040", "MIS-C-M"="#FC9272")

cloned = readChangeoDb("processed/all_cloned_data.tsv")

cloned = filter(cloned, !is.na(c_call))
cloned$major_c = substr(cloned$c_call,1,4)
cloned = filter(cloned, cell_type != "T-NK-B cell doublets")
cloned = filter(cloned, major_c != "IGKC")

cloned[cloned$timepoint == "B" &
	cloned$condition == "COVID19",]$condition = "COVID-B"
cloned[cloned$timepoint == "A" &
	cloned$condition == "COVID19",]$condition = "COVID-A"

cloned = cloned[order(cloned$condition,cloned$sample_id),]

cloned$c_call = factor(cloned$c_call,
	levels = c("IGHM","IGHD","IGHG3","IGHG1","IGHA1",
		"IGHG2","IGHG4","IGHE","IGHA2"))

cloned$major_c = factor(cloned$major_c,
  levels = c("IGHM","IGHD","IGHG","IGHA","IGHE"))

cloned$cell_type = factor(cloned$cell_type)
cloned$bcell_type_combined = cloned$bcell_type
cloned[grepl("lasma",cloned$bcell_type),]$bcell_type_combined = "Plasmablast"
cloned[grepl("emory",cloned$bcell_type),]$bcell_type_combined = "Memory B-cell"
cloned$bcell_type_combined = factor(cloned$bcell_type_combined)
cloned$bcell_type = factor(cloned$bcell_type)

cloned$condition2 = cloned$condition
cloned[cloned$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1") &
  cloned$condition == "MIS-C",]$condition2 = "MIS-C-S"
cloned[!cloned$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1") &
  cloned$condition == "MIS-C",]$condition2 = "MIS-C-M"

plasma_props = cloned %>%
  filter(bcell_type_combined == "Plasmablast") %>%
	group_by(condition,condition2,sample_id,c_call,.drop=FALSE) %>%
	summarize(n=n()) %>%
  mutate(Frequency=n/sum(n),cells=sum(n)) %>%
  ungroup()

if(sum(table(plasma_props$sample_id) != 
    n_distinct(plasma_props$c_call)) != 0){
    stop("grouping failed")
}

boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
# Compare plasmablast frequency across cohorts for different isotypes
my_comparisons <- list( c("MIS-C", "C.HD"))
p1 = ggboxplot(filter(plasma_props,condition %in%
  c("C.HD","MIS-C","MIS-C-R") & c_call == "IGHG1"),
    x = "condition", y = "Frequency", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
	theme(axis.text.x = element_text(angle = 20,
	 vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"
     ))

p2 = ggboxplot(filter(plasma_props,condition %in%
  c("C.HD","MIS-C","MIS-C-R") & c_call == "IGHG3"),
    x = "condition", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"
     ))

ggsave(p1, file = "results/plasma_ig1_frequency_1.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

ggsave(p2, file = "results/plasma_ig3_frequency_1.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

my_comparisons <- list(c("A.HD", "COVID-A"), c("A.HD", "COVID-B"))
p1g1a = ggboxplot(filter(plasma_props,condition %in%
  c("A.HD","COVID-A","COVID-B") & c_call == "IGHG1"),
    x = "condition", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
    theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

p2g3a = ggboxplot(filter(plasma_props,condition %in%
  c("A.HD","COVID-A","COVID-B") & c_call == "IGHG3"),
    x = "condition", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p1g1a, file = "results/plasma_ig1_frequency_adult.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

ggsave(p2g3a, file = "results/plasma_ig3_frequency_adult.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)


# Compare memory B cell proportions across cohorts
memory_props = cloned %>%
  filter(bcell_type_combined == "Memory B-cell") %>%
  group_by(condition,condition2,sample_id,c_call,.drop=FALSE) %>%
  summarize(n=n()) %>%
    mutate(Frequency=n/sum(n),cells=sum(n)) %>%
    ungroup()

my_comparisons <- list( c("MIS-C", "C.HD"))
pmemory = ggboxplot(filter(memory_props,condition %in%
  c("C.HD","MIS-C","MIS-C-R") & c_call == "IGHM"),
    x = "condition", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of memory B cells")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(pmemory, file = "results/plasma_igm_memory_adult.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

# Compare IgG plasmablast frequencies across subtypes
plasma_props_filtered = filter(plasma_props, sample_id != "P6.1")
plasma_props_filtered$condition2 = 
  factor(plasma_props_filtered$condition2,
    levels=c("C.HD","MIS-C-S","MIS-C-M","MIS-C-R"))

my_comparisons <- list( c("MIS-C-S", "C.HD"))
p1g1 = ggboxplot(filter(plasma_props_filtered,condition2 %in%
  c("C.HD","MIS-C-S","MIS-C-M","MIS-C-R") & c_call == "IGHG1"),
    x = "condition2", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols_sm) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

p2g3 = ggboxplot(filter(plasma_props_filtered,condition2 %in%
  c("C.HD","MIS-C-S","MIS-C-M","MIS-C-R") & c_call == "IGHG3"),
    x = "condition2", y = "Frequency",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols_sm) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   facet_wrap(~c_call)+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p1g1, file = "results/plasma_ig1_frequency_ms.pdf", height = 1.5, 
  width =1.2, useDingbats=FALSE)

ggsave(p2g3, file = "results/plasma_ig3_frequency_ms.pdf", height = 1.5, 
  width =1.2, useDingbats=FALSE)

# Plot frequency of unmutated plasma IGHG clones
cell_medians = cloned %>%
	filter(locus == "IGH" & !is.na(c_call)) %>%
	filter(major_c %in% c("IGHG")) %>%
	group_by(condition,condition2,sample_id,bcell_type_combined,clone_id) %>%
	summarize(median_freq = median(mu_freq),
		condition=unique(condition)) %>%
	mutate(mutated = median_freq >= 0.01)

cell_means = cell_medians %>%
	group_by(condition,condition2,sample_id,bcell_type_combined,.drop=FALSE) %>%
	summarize(pmutated = mean(mutated),
		nclones = n())

if(sum(table(cell_means$sample_id) != 
    n_distinct(cell_means$bcell_type_combined)) != 0){
    stop("grouping failed")
}

mutated_plasma = filter(cell_means,
	grepl("Plasma",bcell_type_combined) & !is.na(pmutated))

my_comparisons <- list( c("MIS-C", "C.HD"))
p3 = ggboxplot(filter(mutated_plasma,condition %in%
  c("C.HD","MIS-C","MIS-C-R")),
    x = "condition", y = "pmutated",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Prop. mutated IgG plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p3, file = "results/plasma_mutated_pediatric.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

my_comparisons <- list(c("A.HD", "COVID-A"), c("A.HD", "COVID-B"))
p4 = ggboxplot(filter(mutated_plasma,!condition %in%
  c("C.HD","MIS-C","MIS-C-R")),
    x = "condition", y = "pmutated",
          outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Prop. mutated IgG plasmablasts")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p4, file = "results/plasma_mutated_adult.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

# Plot frequencies of plasmablast isotypes
isotypes = cloned %>%
  filter(bcell_type_combined == "Plasmablast") %>%
  group_by(condition,sample_id,major_c, .drop=FALSE) %>%
  summarize(n = n()) %>%
  mutate(Frequency = n/sum(n))

pisotype = ggplot(filter(isotypes,condition %in% c("C.HD","MIS-C")),
  aes(x=sample_id,y=Frequency,fill=major_c))+
  geom_bar(stat='identity',color="black",size=0.2,width=1)+
  scale_fill_brewer(palette="Dark2",name="Isotype")+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"))+
  theme(legend.key.size = unit(0.5,"line"))+
 theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(pisotype, file = "results/plasma_isotypebar_pediatric.pdf", height = 1.5, 
  width =2.5, useDingbats=FALSE)

aisotype = ggplot(filter(isotypes,!condition %in% c("C.HD","MIS-C","MIS-C-R")),
  aes(x=sample_id,y=Frequency,fill=major_c))+
  geom_bar(stat='identity',color="black",size=0.2,width=1)+
  scale_fill_brewer(palette="Dark2",name="Isotype")+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"))+
  theme(legend.key.size = unit(0.5,"line"))+
 theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(aisotype, file = "results/plasma_isotypebar_adult.pdf", height = 1.5, 
  width =2.5, useDingbats=FALSE)

# Plot subisotype frequencies
subisotypes = cloned %>%
  filter(bcell_type_combined == "Plasmablast") %>%
  group_by(condition,sample_id,c_call, .drop=FALSE) %>%
  summarize(n = n()) %>%
  mutate(Frequency = n/sum(n))

subpisotype = ggplot(filter(subisotypes,condition %in% c("C.HD","MIS-C")),
  aes(x=sample_id,y=Frequency,fill=c_call))+
  geom_bar(stat='identity',color="black",size=0.2,width=1)+
  scale_fill_brewer(palette="Set1",name="Subisotype")+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"))+
  theme(legend.key.size = unit(0.5,"line"))+
  theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(subpisotype, file = "results/plasma_subisotypebar_pediatric.pdf", height = 1.5, 
  width =2.5, useDingbats=FALSE)

subaisotype = ggplot(filter(subisotypes,!condition %in% c("C.HD","MIS-C","MIS-C-R")),
  aes(x=sample_id,y=Frequency,fill=c_call))+
  geom_bar(stat='identity',color="black",size=0.2,width=1)+
  scale_fill_brewer(palette="Set1",name="Subisotype")+
  ylab("Proportion of plasmablasts")+
  xlab("")+
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"))+
  theme(legend.key.size = unit(0.5,"line"))+
  theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(subaisotype, file = "results/plasma_subisotypebar_adult.pdf", height = 1.5, 
  width =2.5, useDingbats=FALSE)

# Clonal diversity analysis of all cells in pediatric cohort
sample_curve_pediatric <- alphaDiversity(
	filter(cloned,condition %in% c("C.HD","MIS-C","MIS-C-R")),
        group="sample_id", clone="clone_id",
		min_q=1.9, max_q=2.1, step_q=0.1,
		ci=0.95, nboot=100, uniform=TRUE)

cq2 = filter(sample_curve_pediatric@diversity, q==2)

cq2$condition = "MIS-C"
cq2[grepl("C.HD",cq2$sample_id),]$condition = "C.HD"
cq2[grepl("\\.2",cq2$sample_id),]$condition = "MIS-C-R"
cq2$condition2 = cq2$condition
cq2[cq2$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1") &
  cq2$condition == "MIS-C",]$condition2 = "MIS-C-S"
cq2[!cq2$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1") &
  cq2$condition == "MIS-C",]$condition2 = "MIS-C-M"

my_comparisons <- list( c("MIS-C", "C.HD"))
p5 = ggboxplot(cq2, x = "condition", y = "d",
          outlier.shape=NA,size=boxsize) +
geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Simpson's diversity")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p5, file = "results/bcell_diversity_pediatric.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

# Clonal diversity analysis of all cells vs covid in adult cohort
sample_curve_adult <- alphaDiversity(
	filter(cloned,condition %in% c("A.HD","COVID-A","COVID-B")),
        group="sample_id", clone="clone_id",
		min_q=1.9, max_q=2.1, step_q=0.1,
		ci=0.95, nboot=100, uniform=TRUE)

aq2 = filter(sample_curve_adult@diversity,q==2)

aq2$condition = "A.HD"
aq2[grepl("\\.1",aq2$sample_id),]$condition = "COVID-A"
aq2[grepl("\\.2",aq2$sample_id),]$condition = "COVID-B"
aq2$condition = factor(aq2$condition,levels=c("A.HD","COVID-A","COVID-B"))
aq2$condition2 = aq2$condition

my_comparisons <- list(c("A.HD", "COVID-B"),c("A.HD", "COVID-A"))
p6 = ggboxplot(aq2, x = "condition", y = "d",
          outlier.shape=NA,size=boxsize) +
geom_jitter(width = 0.05, size= 0.05, height=0, aes(color=condition2))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Simpson's diversity")+
  xlab("")+
  scale_color_manual(values=cols) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
    size = 5, color = "black"), axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20,
   vjust = 0.7, hjust=0.5)) + theme(legend.position = "none") +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
   theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(p6, file = "results/bcell_diversity_adult.pdf", height = 1.5, 
  width =1, useDingbats=FALSE)

cell_props = cloned %>%
  group_by(condition,sample_id,bcell_type_combined,.drop=FALSE) %>%
  summarize(n=n()) %>%
  mutate(Frequency=n/sum(n),cells=sum(n)) %>%
  ungroup()

# Clonal diversity analysis of all cells within MIS-C patients
sample_curve_pediatric2 <- alphaDiversity(
  filter(cloned,condition %in% c("MIS-C")),
        group="sample_id", clone="clone_id",
    min_q=2, max_q=2, step_q=0,
    ci=0.95, nboot=100, uniform=TRUE)

bcr = filter(sample_curve_pediatric2@diversity, q==2)

bcr_range = range(bcr$d)
key = unique(select(cloned,orig.ident,sample_id))
namekey = unlist(key[,2])
names(namekey) = unlist(key[,1])

# B cell subtype frequencies and unmutated frequencies
bcell_medians = cloned %>%
    filter(locus == "IGH" & !is.na(c_call)) %>%
    filter(major_c %in% c("IGHG")) %>%
    group_by(condition,sample_id,bcell_type,clone_id) %>%
    summarize(median_freq = median(mu_freq),
        condition=unique(condition)) %>%
    mutate(mutated = median_freq >= 0.01)

bcell_means = bcell_medians %>%
    group_by(condition,sample_id,bcell_type,.drop=FALSE) %>%
    summarize(pmutated = mean(mutated),
        nclones = n())

bcell_props = cloned %>%
  group_by(condition,sample_id,bcell_type,.drop=FALSE) %>%
  summarize(n=n()) %>%
    mutate(Frequency=n/sum(n),cells=sum(n)) %>%
    ungroup()

# make plot of TCR vs BCR diversity
tcr = read.csv("intermediates/diversity.csv",stringsAsFactors=FALSE)
tcr = filter(tcr,q==2 & group %in% c("MIS-C"))

tcr = filter(tcr,sample %in% names(namekey))
tcr$sample_id = namekey[tcr$sample]

tcr$tcr_d = tcr$d
m = match(tcr$sample_id,bcr$sample_id)

tcr$bcr_d = bcr[m,]$d
tcr$bcr_sample_id = bcr[m,]$sample_id
mean(tcr$bcr_sample_id == tcr$sample_id)

tcr$Condition = tcr$group
tcr[tcr$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1"),]$Condition = "MIS-C-S"
tcr[!tcr$sample_id %in% c("P1.1", "P2.1", "P3.1", "P7.1"),]$Condition = "MIS-C-M"
cols = c("MIS-C-S"="#C94040", "MIS-C-M"="#FC9272")
shapes = c("MIS-C-S"=21, "MIS-C-M"=24)

tcr = filter(tcr, sample_id != "P6.1")

bcr_tcr = ggplot(filter(tcr,condition=="ki67+memory_CD4"),
    aes(y=tcr_d, x=bcr_d, fill=Condition))+
geom_point(aes(shape=Condition),stroke=0.25)+facet_wrap(~condition,nrow=2,
    scales="free",labeller = labeller(condition = c("ki67+memory_CD4"="KI67+ Memory CD4")))+
scale_shape_manual(values=shapes)+
theme_bw()+scale_fill_manual(values=cols)+
xlab("Total B cell diversity")+ylab("T cell diversity")+
theme(strip.background =element_rect(fill="white"))+
theme_classic(base_size = 5) +
theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
  size = 5, color = "black"), axis.text.y = element_text(color = "black")) +
theme(axis.text.x = element_text(angle = 20,vjust = 0.7, hjust=0.5)) + 
 scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
 theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
 xlim(bcr_range[1],bcr_range[2])+
 theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(bcr_tcr,file="results/bcr_tcr.pdf",width=1.75,
  height=1.5,useDingbats=FALSE)

# make frequency, % unmutated, and diversity correlation plots
comp1_plasma = filter(cell_means,(condition %in% c("MIS-C")) &
    grepl("lasma",bcell_type_combined) & !is.na(pmutated))
comp2_plasma = filter(cell_props,(condition %in% c("MIS-C")) &
    grepl("lasma",bcell_type_combined) & !is.na(Frequency))

m = match(comp1_plasma$sample_id, bcr$sample_id)
sum(comp1_plasma$sample_id != bcr[m,]$sample_id)
comp1_plasma$diversity = bcr[m,]$d

m = match(comp1_plasma$sample_id, comp2_plasma$sample_id)
sum(comp1_plasma$sample_id != comp2_plasma[m,]$sample_id)
comp1_plasma$Frequency = comp2_plasma[m,]$Frequency

comp1_plasma$Condition = "MIS-C-M"
comp1_plasma[comp1_plasma$sample_id %in% 
    c("P1.1", "P2.1", "P3.1", "P7.1"),]$Condition = "MIS-C-S"

everything = ggplot(comp1_plasma,aes(x=diversity,y=Frequency,
    fill=pmutated))+
geom_point(aes(shape=Condition),stroke=0.25)+
theme_bw()+
scale_shape_manual(values=shapes)+
xlab("Total B cell diversity")+ylab("Plasmablast frequency")+
theme(strip.background =element_rect(fill="white"))+
scale_fill_viridis_c(option="inferno",direction=1,
  name="Prop. Mutated\nIgG clones")+
theme_classic(base_size = 5) +
theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, 
  size = 5, color = "black"), axis.text.y = element_text(color = "black")) +
theme(axis.text.x = element_text(angle = 20,vjust = 0.7, hjust=0.5)) + 
 scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
 theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
 xlim(bcr_range[1],bcr_range[2])+
 theme(axis.line = element_line(colour = 'black', size = axissize),
    strip.background = element_rect(
     color="black", fill="white", size=panelsize, linetype="solid"))

ggsave(everything,file="results/diversity_v_plasmablast_all.pdf",width=1.75,
  height=1.5,useDingbats=FALSE)
