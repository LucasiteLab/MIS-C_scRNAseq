library(ggplot2)
library(dplyr)
library(ggrepel)


volcano <- read.csv("Sheets/all_de_proteins_with_myeloid_and_liver_status.csv")
volcano[,1] <- NULL


sol <- c("IL1RN", "IL1R1", "SAA1", "IL1RL1", "HAMP", "CXCL10", "CHI3L1", # annotated soluble genes
         "TNFRSF1B", "TNFRSF1B", "PTGDS", "IL1R2", "TCN2", "CXCL9", "ITIH1", "FCGR3B", 
         "INHBA", "CCDC126", "CD163", "CXCL11", "FAM3B", "CCL23", "LILRA5", "C1QC", 
         "SCGB3A1", "INHBB", "ADAMTS13", "ABO", "CCL23", "HGF", "BMP4", "APOC1", "APOA1", 
         "INHBA", "QPCT", "RNASE6", "KNG1", "SERPINA1", "NPTX2", "GSN", "CXCL6", "PLTP", 
         "HBEGF", "RBP4", "TFPI", "PTK7", "HP", "SERPINC1", "VEGFA", "TGFB1") 

pm <-c("MRC1", "IL1R1", "IL1RL1", "GPNMB", "CDON", "GPNMB", "TNFRSF1B",  # annotated plasma membrane associated
       "IL1R2", "FCGR3B", "CD163", "RET", "LILRA5", "ROR1", "FLRT2", "SIGLEC5", "IL2RA",
       "LILRB2", "OMG", "LCT", "PLXNC1", "UNC5B", "CD300C", "MMP17", "SIGLEC14", "HAVCR2", 
       "TNFRSF8", "CLEC7A", "SERPINA1", "CD177", "UNC5B", "FCGR1A", "RSPO2", "HBEGF",
       "CLEC12A", "PTK7", "SERPINC1", "PECAM1", "ICAM3", "CLEC1B", "FCGR2B", 
       "NTM", "ESAM", "SCARA5") 

soluble <- intersect(sol, pm)

volcano_soluble <- volcano %>% filter(EntrezGeneSymbol %in% soluble)

volcano_soluble <- volcano_soluble %>% filter(!(AptName %in% c("TNFRSF1B.3152.57", "LILRA5.8766.29")))

soluble <- volcano_soluble$EntrezGeneSymbol
soluble <- c(soluble, "SELE")

volcano$labeling <- rep(NA, nrow(volcano))

for(i in 1:nrow(volcano)){
  if(volcano[i,14] == "True" & abs(volcano[i,4]) > 1 & -log10(volcano[i,8]) > 1.3010){
    volcano[i,16] <- "Myeloid-derived soluble factor"
  } else if(volcano[i,13] == "True" & abs(volcano[i,4]) > 1 & -log10(volcano[i,8]) > 1.3010) {
    volcano[i,16] <- "Myeloid-derived plasma membrane"
  } else if(volcano[i,3] == "SELE") {
    volcano[i,16] <- "E-selectin"
  } else if(abs(volcano[i,4]) > 1 & volcano[i,12] > -log10(0.05)) {
    volcano[i,16] <- "Pass"
  } else {
    volcano[i,16] <- "All other"
  }
}



plot1 <- ggplot(volcano, aes(x = logFC, y = log.adj.P.Val)) + 
  geom_point(data = subset(volcano, labeling == 'All other'), aes(color = labeling), size =3) +
  geom_point(data = subset(volcano, labeling == 'Pass'), aes(color = labeling), size =3) +
  geom_point(data = subset(volcano, labeling == 'Myeloid-derived plasma membrane'), aes(color = labeling), size =3) +
  geom_point(data = subset(volcano, labeling == 'Myeloid-derived soluble factor'), aes(color = labeling), size =3) +
  geom_point(data = subset(volcano, labeling == 'E-selectin'), aes(color = labeling), size = 4) +
  geom_text_repel(aes(label=ifelse(EntrezGeneSymbol %in% soluble, EntrezGeneSymbol ,'')),
                  box.padding =0.75) +
  geom_text_repel(aes(label=ifelse(logFC > 3.5, EntrezGeneSymbol ,'')),
                  box.padding =0.75) +
  scale_color_manual(values = c("#d9d9d9", "#cb181d", "#35B779FF", "#F7D03CFF", "#1f1d75")) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  ylab("-log10(p.adj)") +
  xlab("logFC") +
  theme(legend.text = element_text(size = 12), 
        axis.line = element_line(colour = 'black', size = 0.6),
        axis.text = element_text(colour = "black")) 
plot1

ggsave(plot1, file = "myeloid_volcano_updated.pdf", height =7 , width = 9)



significant <- volcano %>% filter(abs(logFC) > 1 & log.adj.P.Val > -log10(0.05))
ns <- volcano %>% filter(abs(logFC) <= 1 | log.adj.P.Val <= -log10(0.05))

# myeloid
myeloid <- volcano %>% filter(myeloid_plasma == "True" | myeloid_soluble == "True")
myeloid_not <- volcano %>% filter(!(myeloid_plasma == "True" | myeloid_soluble == "True"))


# significant 
sig_genes <- unique(significant$AptName)
# non-significant  
ns_genes <- unique(ns$AptName)
# myeloid 
myeloid_genes <- unique(myeloid$AptName)
# not myeloid 
nm_genes <- unique(myeloid_not$AptName)
# myeloid AND DE 42
myeloid_and_de <- intersect(sig_genes, myeloid_genes)
# myeloid and not DE 247
myeloid_and_nde <- intersect(ns_genes, myeloid_genes)
# not myeloid and DE 166
nm_and_de <- intersect(sig_genes, nm_genes)
# not myeloid and not DE 4251
nde_and_nm <- intersect(ns_genes, nm_genes)

ct <- data.frame("DE" = c(2,206), "NDE" = c(8,4490)) # values change depending on gene set of interest
rownames(ct) <- c("myeloid", "not_myeloid")
ct <- as.matrix(ct)

fisher <- fisher.test(ct, alternative = "greater")
