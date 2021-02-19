library(ggplot2)
library(tidyverse)
library(ggsignif)

cd4_cd8 <- read.csv("Sheets/GeneUsage_SampleLevel.csv")
cd4_cd8[,1] <- NULL

cd4_cd8_vb11.2 <- cd4_cd8 %>% filter(Genes == "TRBV11-2")

cd4_vb11.2 <- cd4_cd8_vb11.2 %>% filter(cell_type == "CD4_Memory")

cd8_vb11.2 <- cd4_cd8_vb11.2 %>% filter(cell_type == "CD8_Memory")

level_order <- c("C.HD", "MIS-C-S", "MIS-C-M")
cd4_vb11.2$Group <- factor(cd4_vb11.2$Group, level = level_order)
cd8_vb11.2$Group <- factor(cd8_vb11.2$Group, level = level_order)


## CD4 memory
# Sig

chd_tmp <- cd4_vb11.2 %>% filter(Group == "C.HD")
miscm_tmp <- cd4_vb11.2 %>% filter(Group == "MIS-C-M")
miscs_tmp <- cd4_vb11.2 %>% filter(Group == "MIS-C-S")


# Example

w.test <- wilcox.test(x = chd_tmp$Freq, y = miscm_tmp$Freq, alternative =  "less", correct = FALSE, exact = FALSE)
pval_mod_cd4 <- w.test$p.value

w.test <- wilcox.test(x = chd_tmp$Freq, y = miscs_tmp$Freq, alternative =  "less", correct = FALSE, exact = FALSE)
pval_severe_cd4 <- w.test$p.value


cols <- c("#6baed6", "#c94040", "#FC9272")
# Plotting

plot1 <- ggplot(cd4_vb11.2, aes(x = Group, y = Freq)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(aes(colour = factor(Group, level = level_order)), size = 0.5, width = 0.12) +
  ggtitle("CD4") +
  xlab("") +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C-S")), annotation = "0.04", 
              size = 0.12, textsize = 2,
              y_position =0.25)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylim(c(0, 0.28)) +
  ylab("TRBV11-2 Frequency")

plot1

ggsave(plot1, file = "CD4_memory_TRBV11-2.pdf", height = 1.5, width =0.8)


### CD8

# Sig

chd_tmp_cd8 <- cd8_vb11.2 %>% filter(Group == "C.HD")
miscm_tmp_cd8 <- cd8_vb11.2 %>% filter(Group == "MIS-C-M")
miscs_tmp_cd8 <- cd8_vb11.2 %>% filter(Group == "MIS-C-S")


# Example

w.test <- wilcox.test(x = chd_tmp_cd8$Freq, y = miscm_tmp_cd8$Freq, alternative = "less", correct = FALSE, exact = FALSE)
pval_mod_cd8 <- w.test$p.value

w.test <- wilcox.test(x = chd_tmp_cd8$Freq, y = miscs_tmp_cd8$Freq, alternative = "less", correct = FALSE, exact = FALSE)
pval_severe_cd8 <- w.test$p.value


# Plotting

plot1 <- ggplot(cd8_vb11.2, aes(x = Group, y = Freq)) +
  geom_boxplot(lwd=0.15, outlier.shape = NA) +
  geom_jitter(aes(colour = factor(Group, level = level_order)), size = 0.5, width = 0.12) +
  ggtitle("CD8") +
  xlab("") +
  scale_color_manual(values = cols) +
  geom_signif(comparisons = list(c("C.HD", "MIS-C-S")), annotation = "0.005", 
              size = 0.12, textsize = 2,
              y_position =0.67)+
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 6, color = "black"),
        axis.text.y = element_text(color = "black"), axis.line = element_line(size = 0.15), 
        legend.position = "none", axis.ticks = element_line(size = 0.15)) +
  ylim(c(0, 0.7)) +
  ylab("TRBV11-2 Frequency")

plot1

ggsave(plot1, file = "CD8_memory_TRBV11-2.pdf", height = 1.5, width =0.8)