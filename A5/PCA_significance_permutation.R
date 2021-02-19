# Kenneth B. Hoehn
# 2/2/2021
# Assess whether cells within conditions are more closely
# clustered in PCA plots than expected by chance

library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(scales)


# Compuate mean intra and inter-cluster
# Euclidian distances using PC1 and PC2
pca_dist = function(x){
  dists = tibble()
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
        x1 = x[i,]$PC1
        y1 = x[i,]$PC2
        x2 = x[j,]$PC1
        y2 = x[j,]$PC2
        clustered = "within"
        if(x[i,]$Condition != x[j,]$Condition){
          clustered = "between"
        }
        distance = sqrt((x1-x2)^2 + (y1-y2)^2)
        dists = bind_rows(dists,tibble(from=x[i,]$sample_id,to=x[j,]$sample_id,
          distance=distance,clustered=clustered))
    }
  }
  dists %>% 
        group_by(clustered) %>% 
        summarize(distance=mean(distance)) %>% 
        spread(clustered,distance) %>%
        mutate(ratio = within/between)
}

reps = 1000

# TCR PCAs
tcr_naive = read.csv("intermediates/PCA_matrix_CD4_Naive.csv")
tcr_naive$celltype = "CD4_Naive"
tcr_memory = read.csv("intermediates/PCA_matrix_CD4_Memory+Ki67.csv")
tcr_memory$celltype = "CD4_Memory+Ki67"

tcr = bind_rows(tcr_naive, tcr_memory)
tcr = select(tcr,-X)
tcr = rename(tcr,sample_id=Samples,Condition=Group)

# BCR PCAs
bcr = read.csv("intermediates/PCA_matrix_bcrs.csv")
bcr = select(bcr,Condition,sample_id,celltype,PC1,PC2,-X)

total = bind_rows(bcr,tcr)

cluster_plots = list()
for(celltype in unique(total$celltype)){
    print(celltype)
    temp = filter(total,!!celltype==celltype)

    # Calculate intra vs inter cluster distances
    dists = pca_dist(temp)
    dists$type = "observed"

    # Calculate intra vs inter cluster distances over
    # the specified number of random permutations
    for(i in 1:reps){
      temp$Condition = sample(temp$Condition)
      perm = pca_dist(temp)
      perm$type = "permuted"
      dists = bind_rows(dists,perm)
    }

    # P value is proportion of permuted repetitions with an intra/inter
    # cluster distance at least as small as observed
    pvalue = signif(mean(filter(dists, type=="permuted")$ratio <= 
        filter(dists, type=="observed")$ratio), digits=2)

    cluster_plots[[celltype]] = ggplot(filter(dists,type=="permuted"),
      aes(x=ratio))+geom_histogram()+
    geom_vline(xintercept=filter(dists,type=="observed")$ratio,color="red")+
    theme_bw()+ggtitle(paste(celltype,", p value:",pvalue))+
    xlab("Mean intra/inter-cluster distance")+
    ylab("Number of repetitions")

}

pdf("results/PCA_significance.pdf",width=8,height=10)
grid.arrange(grobs=cluster_plots,nrow=3)
dev.off()