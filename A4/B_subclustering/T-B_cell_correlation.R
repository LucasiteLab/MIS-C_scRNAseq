# cell correlation between Ki67+CD4+Tcells and PBs
cells <- read.csv("cell_proportion_10.txt", header=T, sep="\t") # cell proportions for scRNAseq

p <- ggplot(aes(x=Ki67CD4Tcells, y=Ki67PBs, color=gps), data = cells) + geom_point(size = 3.5)
p <- p + geom_smooth(method = "lm", color="black", se=TRUE)
p <- p + theme_classic()
p <- p + scale_color_manual(name="gps",values=c("MIS-C"="#c94040"))
p

summary(lm(cells[,2] ~ cells[,1]))


# Flow data

library(ggplot2)
library(ggrepel)

cells <- read.csv("Correlation_flow.csv")
cells <- cells[1:14,]
rownames(cells) <- cells[,1]
cells[,1] <- NULL

summary(lm(cells[,2] ~ cells[,1]))
cor.test(cells[,1], cells[,2], method = "pearson") 

names(cells) <- c("tcell", "bcell", "names", "Severity")

cols <- c("#c94040", "#FC9272")
cells$Severity <- factor(cells$Severity, level = c("MIS-C-S", "MIS-C-M"))


p <- ggplot(aes(x=tcell, y=bcell, color=Severity), data = cells) +
      geom_smooth(method = "lm", color="gray30", se=TRUE, fill = "gray70", size = 0.2) + 
      geom_point(size = 0.5) + theme_classic(base_size = 7) +
      geom_text_repel(aes(label = names), color = "black", size = 1, vjust = 0) +
      scale_color_manual(values = cols) +
      ylab("%CD19+CD27++CD38++KI67+/CD19+") +
      xlab("%CD4+KI67+/CD3+") +
      theme(axis.line = element_line(size = 0.15), legend.position = "none", 
            axis.text = element_text(color = "black"))
p

ggsave(p, file = "correlation_flow_wlabels.pdf", height = 2, width = 2.2)


