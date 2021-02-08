library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(uwot)
library(vizier)
library(MASS)
library(grid)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
p_c<-list()
g_c<-list()


data.path<-"."
#Input the information for "condition","sample"
anno.files<-c("Annotation_Vgene_10x_CD4_Memory.csv","Annotation_Vgene_10x_CD8_Memory.csv")
#Input the corresponding indicator matrix for gene-usage
count.files<-c("Count_Vgene_10x_CD4_Memory.csv","Count_Vgene_10x_CD8_Memory.csv")

for (i in 1:length(anno.files)) {
  set.seed(100)
  anno.df<-read.csv(paste0(data.path,anno.files[i]),header = T)
  anno.df$gene_type<-substr(anno.files[i],12,12)
  anno.df$cell_type<-substr(anno.files[i],22,(nchar(anno.files[i])-4))

  
  count.df<-read.csv(paste0(data.path,count.files[i]),header = T)
  count.df$cell_id<-rownames(count.df)
  tmp.df<-merge(anno.df,count.df,by="cell_id")
  tmp.df$condition[which(tmp.df$condition=="Mild MIS-C")]="MIS-C-M"
  tmp.df$condition[which(tmp.df$condition=="Severe MIS-C")]="MIS-C-S"
  rownames(tmp.df)<-tmp.df$cell_id
  full.df<-tmp.df[-which(tmp.df$sample=="Y127-1"),] #Filter out bad-quality sample

  #Sample-level gene-usage frequency
  sample.df<-full.df[,-which(names(full.df)%in%c("cell_id"))] %>%
    group_by(condition,sample, gene_type,cell_type) %>%
    summarise_each(funs(sum))

  #For selected groups of patients in the data
  group_comb<-list(c("C.HD","MIS-C-S","MIS-C-M"),c("A.HD","COVID19-A"),c("A.HD","COVID19-B"))
  pca_group_figure<-list()
  for(g in 1:length(group_comb)){
    select.df<-sample.df[which(sample.df$condition%in%group_comb[[g]]),]
    select.df_clean<-cbind(select.df[,1:4],select.df[,-c(1:4)][,colSums(select.df[,-c(1:4)])!=0])
    sample_fraction.df<-select.df_clean
    sample_fraction.df[,-c(1:4)]<-sample_fraction.df[,-c(1:4)]/matrix(rep(rowSums(sample_fraction.df[,-c(1:4)]),ncol(sample_fraction.df)),ncol=ncol(sample_fraction.df))
    colnames(sample_fraction.df)<-gsub("[.]", "-",  colnames(sample_fraction.df))
    pca_frac<-prcomp(sample_fraction.df[,-c(1:4)],scale. = T,center = T)

    #Figure: PCA for different groups of patients
    tmp_pca_figure<-ggplot(cbind(as.data.frame(sample_fraction.df),pca_frac$x),aes(x=PC1, y=PC2,color=condition))+
      scale_color_manual(values=c("MIS-C-M"="#FC9272","C.HD"="#6BAED6", "MIS-C-S"="#C94040", "MIS-C-R"="#969696", "A.HD"="#9970AB", "COVID19-A"="#EC7014", "COVID19-B"="#FEC44F") )+theme_classic()+
      theme(plot.title = element_text(hjust = 0.5,size=7),legend.position  = "none",axis.title.x = element_text(size=6),axis.title.y = element_text(size=6),axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))+xlab(paste0("PC1 (",percent(summary(pca_frac)[["importance"]]["Proportion of Variance","PC1"],accuracy = 0.01),")"))+
      ylab(paste0("PC2 (",percent(summary(pca_frac)[["importance"]]["Proportion of Variance","PC2"],accuracy = 0.01),")"))+geom_point(size=0.2)+ggtitle("PCA TRBV-Usage")
    q <- ggplot_build(tmp_pca_figure)
    pca_group_figure[[g]] <- ggplot_gtable(q)

    
    if(identical(group_comb[[g]],c("C.HD","MIS-C-S","MIS-C-M"))){
      sample.df_long<-data.frame(Samples=factor(rep(sample_fraction.df$sample,(ncol(sample_fraction.df)-4))),Group=factor(rep(sample_fraction.df$condition,(ncol(sample_fraction.df)-4))),gene_type=rep(sample_fraction.df$gene_type,(ncol(sample_fraction.df)-4)),cell_type=rep(sample_fraction.df$cell_type,(ncol(sample_fraction.df)-4)),Genes=factor(rep(colnames(sample_fraction.df)[-c(1:4)],each=nrow(sample_fraction.df))),Freq=as.vector(as.matrix(sample_fraction.df[,-c(1:4)])))
      group.df <- data_summary(sample.df_long, varname="Freq", 
                               groupnames=c("Group", "Genes"))
      reorder_levels_2=group.df$Genes[which(group.df$Group=="MIS-C-S")][order(group.df$Freq[which(group.df$Group=="MIS-C-S")]-group.df$Freq[which(group.df$Group=="C.HD")],decreasing = F)]
      
      group.df$Genes <- factor(group.df$Genes, levels = reorder_levels_2)
      
      #Figure: Gene-Usage differences
      tmp_usage <- ggplot(group.df, aes(x=Genes, y=Freq, fill=Group)) + 
        geom_bar(stat="identity", position=position_dodge()) +ylab("Frequency")+
        geom_errorbar(aes(ymin=Freq, ymax=Freq+sd), width=.2,
                      position=position_dodge(.9))+
        scale_fill_manual(values=c("MIS-C-M"="#FC9272","C.HD"="#6BAED6", "MIS-C-S"="#C94040", "MIS-C-R"="#969696", "A.HD"="#9970AB", "COVID19-A"="#EC7014", "COVID19-B"="#FEC44F") )+
        coord_flip()+ 
        theme(plot.title = element_blank(),axis.title.x = element_text(size=6),legend.position ="none" ,axis.title.y = element_blank(),axis.text.x = element_text(size=6),axis.text.y = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      p <- ggplot_build(tmp_usage)
      g_c[[i]] <- ggplot_gtable(p)
    }
  }
  p_c[[i]]<-pca_group_figure
}

for(ct in 1:2){
  #Save the PCA plots
  if(ct==1){celltype="CD4_Memory"}else if(ct==2){celltype="CD8_Memory"}
  pdf(paste0(celltype,"_PCA.pdf"),width = 1.5,height = 1.5)
  for(i in 1:length(p_c[[ct]])){
    grid.draw(p_c[[ct]][[i]])
    if(i!=length(p_c[[ct]])){grid.newpage()}
  }
  dev.off()
  
  #Save the gene-usage plot
  pdf(paste0(celltype,"_GeneUsage.pdf"),width = 4,height = 6)
  grid.draw(g_c[[ct]])
  dev.off()
}


