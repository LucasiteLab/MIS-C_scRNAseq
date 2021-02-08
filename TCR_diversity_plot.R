# generate figure S4.i

library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(1)

home.dir = 'C:\\users\\hill103\\Documents\\GitHub\\Yale-MISC-PBMC\\RawData'
group.mapping = c( 'Y111-1' =	'MIS-C',
                   'Y113-1' = 'MIS-C',
                   'HD_32M' = 'A.HD',
                   'HD_35F' =	'A.HD',
                   'HD_36M' = 'A.HD',
                   'NS0A'   =	'COVID19-A',
                   'NS0B'   =	'COVID19-B',
                   'NS1A'   =	'COVID19-A',
                   'NS1B'   =	'COVID19-B',
                   'TS2A'   =	'COVID19-A',
                   'TS2B'   =	'COVID19-B',
                   'TS3A'   = 'COVID19-A',
                   'TS3B'   =	'COVID19-B',
                   'TP9B'   =	'COVID19-B',
                   'TP8B'	  = 'COVID19-B',
                   'C27'    =	'A.HD',
                   'C32'    =	'A.HD',
                   'C33'    =	'A.HD',
                   'C39'    =	'A.HD',
                   'HA5876' = 'A.HD',
                   'HA5877'	= 'A.HD',
                   'HA5894' =	'A.HD',
                   'HA5952'	= 'A.HD',
                   'HA5953' = 'A.HD',
                   'HA5957' = 'A.HD',
                   'Y117-1' = 'MIS-C',
                   'Y124-1' = 'MIS-C',
                   'Y125-1' =	'MIS-C',
                   'Y127-1' =	'MIS-C',
                   'Y129-1' =	'MIS-C',
                   'Y117-R' = 'MIS-C-R',
                   'Y124-R'	= 'MIS-C-R',
                   'NC-13F' = 'C.HD',
                   'Y54-4'  = 'C.HD',
                   'Y70-4'  = 'C.HD',
                   'Y28-2'  =	'C.HD',
                   'Y28-4'  =	'C.HD',
                   'Y29-2'  =	'C.HD'
)
colors = c('A.HD'      = '#9970AB',
           'C.HD'      = '#6BAED6',
           'COVID19-A' = '#EC7014',
           'COVID19-B' = '#FEC44F',
           'MIS-C'     = '#C94040',
           'MIS-C-R'   = '#969696')
misc.mapping =  c('Y111-1' = 'Severe MIS-C',
                  'Y113-1' = 'Severe MIS-C',
                  'Y117-1' = 'Severe MIS-C',
                  'Y129-1' = 'Severe MIS-C',
                  'Y124-1' = 'Mild MIS-C',
                  'Y125-1' = 'Mild MIS-C',
                  'Y127-1' = 'Mild MIS-C',
                  'NC-13F' = 'C.HD',
                  'Y54-4'  = 'C.HD',
                  'Y70-4'  = 'C.HD',
                  'Y28-2'  = 'C.HD',
                  'Y28-4'  = 'C.HD',
                  'Y29-2'  = 'C.HD',
                  'HD_32M' = 'A.HD',
                  'HD_35F' = 'A.HD',
                  'HD_36M' = 'A.HD',
                  'NS0A'   = 'COVID19-A',
                  'NS0B'   = 'COVID19-B',
                  'NS1A'   = 'COVID19-A',
                  'NS1B'   = 'COVID19-B',
                  'TS2A'   = 'COVID19-A',
                  'TS2B'   = 'COVID19-B',
                  'TS3A'   = 'COVID19-A',
                  'TS3B'   = 'COVID19-B',
                  'TP9B'   = 'COVID19-B',
                  'TP8B'	 = 'COVID19-B',
                  'C27'    = 'A.HD',
                  'C32'    = 'A.HD',
                  'C33'    = 'A.HD',
                  'C39'    = 'A.HD',
                  'HA5876' = 'A.HD',
                  'HA5877' = 'A.HD',
                  'HA5894' = 'A.HD',
                  'HA5952' = 'A.HD',
                  'HA5953' = 'A.HD',
                  'HA5957' = 'A.HD',
                  'Y117-R' = 'MIS-C-R',
                  'Y124-R' = 'MIS-C-R')
misc.colors = c('A.HD'         = '#9970AB',
                'C.HD'         = '#6BAED6',
                'COVID19-A'    = '#EC7014',
                'COVID19-B'    = '#FEC44F',
                'Severe MIS-C' = '#C94040',
                'Mild MIS-C'   = '#FC9272',
                'MIS-C-R'      = '#969696')


# Read diversity data

diversity = read.table(file = file.path(home.dir, 'diversity.csv'), header = TRUE, sep = ',', stringsAsFactors = F)

# Update COVID19 group
for (i in 1:nrow(diversity)) {
  if (diversity[i, 'group'] == 'COVID19') {
    # The last character of sample name is A or B
    this.sample = diversity[i, 'sample']
    if (substr(this.sample, nchar(this.sample), nchar(this.sample)) == 'A') {
      diversity[i, 'group'] = 'COVID19-A'
    }
    else if (substr(this.sample, nchar(this.sample), nchar(this.sample)) == 'B') {
      diversity[i, 'group'] = 'COVID19-B'
    }
    else {stop('Time points is not A nor B!')}
  }
}

# add severe / mild MISC annotation
diversity$misc_group = misc.mapping[diversity$sample]


# Define function

oneDiversityPlot = function(diversity.matrix, comps) {
  # plot diversity boxplot
  
  index.names = c("Richness", "Shannon", "Simpson", "Shannon/Richness", "Simpson/Richness")
  p.list = list()
  i = 1
  p.list[[length(p.list)+1]] = diversity.matrix[diversity.matrix$q==i-1, ] %>%
    group_by(group) %>%
    #mutate(outlier = ifelse(is_outlier(d), sample, NA)) %>%
    ggplot(., aes(x = group, y = d, label = sample)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.25) +
    geom_jitter(mapping=aes(color=misc_group), position = position_jitter(seed = 1), size = 1) +
    labs(x = '', y = index.names[i]) +
    scale_color_manual(values = misc.colors) +
    scale_x_discrete() +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text = element_text(size=10, family="Helvetica", color="black"),
          axis.text.x = element_text(angle=45, hjust=1, vjust=0.9),
          axis.title = element_text(size=10),
          axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35)) +
    stat_compare_means(comparisons=comps, method="wilcox.test", label="p.format", size=3)
  #geom_text(position = position_jitter(seed = 1), hjust = -0.1, size = 3) + 
  
  i = 4
  temp.matrix = diversity.matrix[diversity.matrix$q==i-3, ]
  temp.matrix.1 = diversity.matrix[diversity.matrix$q==0, ]
  temp.matrix[, 'd'] = temp.matrix[, 'd'] / temp.matrix.1[, 'd']
  p.list[[length(p.list)+1]] = temp.matrix %>%
    group_by(group) %>%
    #mutate(outlier = ifelse(is_outlier(d), sample, NA)) %>%
    ggplot(., aes(x = group, y = d, label = sample)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.25) +
    geom_jitter(mapping=aes(color=misc_group), position = position_jitter(seed = 1), size = 1) +
    labs(x = '', y = index.names[i]) +
    scale_color_manual(values = misc.colors) +
    scale_x_discrete() +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text = element_text(size=10,  family="Helvetica", color="black"),
          axis.text.x = element_text(angle=45, hjust=1, vjust=0.9),
          axis.title = element_text(size=10),
          axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35)) +
    stat_compare_means(comparisons=comps, method="wilcox.test", label="p.format", size=3)
  #geom_text(position = position_jitter(seed = 1), hjust = -0.1, size = 3)
  
  return(p.list)
}

# Figure C: diversity analysis plot (Ki67 and Memory T cells)

q.list = list()


## Group C.HD, MIS-C, MIS-C-R

### CD4+ T cells (excluding naive)


tmp.diversity = diversity[(diversity$condition=='ki67+memory_CD4') & (diversity$group %in% c('C.HD', 'MIS-C', 'MIS-C-R')), ]

q.list = c(q.list, oneDiversityPlot(tmp.diversity, list(c('C.HD', 'MIS-C'))))


### CD8+ T cells (excluding naive)


tmp.diversity = diversity[(diversity$condition=='ki67+memory_CD8') & (diversity$group %in% c('C.HD', 'MIS-C', 'MIS-C-R')), ]

q.list = c(q.list, oneDiversityPlot(tmp.diversity, list(c('C.HD', 'MIS-C'))))


## Group A.HD, COVID19-A, COVID19-B

### CD4+ T cells (excluding naive)


tmp.diversity = diversity[(diversity$condition=='ki67+memory_CD4') & (diversity$group %in% c('A.HD', 'COVID19-A', 'COVID19-B')), ]

q.list = c(q.list, oneDiversityPlot(tmp.diversity, list(c('A.HD', 'COVID19-A'), c('A.HD', 'COVID19-B'))))


### CD8+ T cells (excluding naive)


tmp.diversity = diversity[(diversity$condition=='ki67+memory_CD8') & (diversity$group %in% c('A.HD', 'COVID19-A', 'COVID19-B')), ]

q.list = c(q.list, oneDiversityPlot(tmp.diversity, list(c('A.HD', 'COVID19-A'), c('A.HD', 'COVID19-B'))))


## Generate plot

#output.filepath = file.path(home.dir, "fig_diversity_plot_ki67+memory.eps")
#setEPS()
# (200 * 200 pt) * 2
#postscript(output.filepath, width=11.12, height=3.95)

# add a extra legend to show all the groups
ggarrange(plotlist=q.list, ncol=8, nrow=1, common.legend = TRUE, legend="bottom",
          legend.grob=cowplot::get_legend(
            diversity[diversity$q==0, ] %>%
              group_by(group) %>%
              ggplot(., aes(x = group, y = d, label = sample)) +
              geom_boxplot(outlier.shape = NA, lwd = 0.25) +
              geom_jitter(mapping=aes(color=misc_group), position = position_jitter(seed = 1)) +
              scale_color_manual(values = misc.colors, breaks = c("C.HD", "Severe MIS-C", "Mild MIS-C", "MIS-C-R", "A.HD", "COVID19-A", "COVID19-B")) +
              theme(legend.position = "bottom",
                    legend.title = element_blank(),
                    legend.text = element_text(size=10),
                    legend.key.size = unit(1.2, "line")) +
              guides(color = guide_legend(nrow = 1))
          ))

#dev.off()

