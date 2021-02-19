
library(dplyr)
library(ggplot2)
library(alakazam)
library(gridExtra)

set.seed(1)

home.dir = '/home/ningshan/Documents/MISC'
project.folders = c('lucas_final', 'ha_final', 'c_final', 'covid_final')
group.mapping = c( 'Y111-1' =	'MIS-C',
                   'Y113-1' = 'MIS-C',
                   'HD_32M' = 'A.HD',
                   'HD_35F' =	'A.HD',
                   'HD_36M' = 'A.HD',
                   'NS0A'   =	'COVID19',
                   'NS0B'   =	'COVID19',
                   'NS1A'   =	'COVID19',
                   'NS1B'   =	'COVID19',
                   'TS2A'   =	'COVID19',
                   'TS2B'   =	'COVID19',
                   'TS3A'   = 'COVID19',
                   'TS3B'   =	'COVID19',
                   'TP9B'   =	'COVID19',
                   'TP8B'	  = 'COVID19',
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
colors = c('A.HD'    = '#9970AB',
           'C.HD'    = '#6BAED6',
           'COVID19' = '#EC7014',
           'MIS-C'   = '#C94040',
           'MIS-C-R' = '#FC9272')


# Loading data
* Use **beta chains** only


meta.data = data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("sample", "group", "color"))))
raw.data = list()

for (folder in project.folders) {
  parent.folder = file.path(home.dir, folder)
  
  child.folders = list.dirs(path = parent.folder, full.names = F, recursive = F)
  
  print(sprintf('Total %d samples under folder %s', length(child.folders), parent.folder))
  
  # iterate all subfolders
  for (f in child.folders) {
    raw.data[[length(raw.data)+1]] = read.table(file = file.path(parent.folder, f, 'filtered_contig_combined_productive-T.tsv'), header = TRUE, sep = '\t', stringsAsFactors = F) %>%
      filter(locus == "TRB") %>%
      # replace cell type '' with 'Missing'
      mutate(annotation = replace(annotation, annotation=='', 'Missing'))
    this.group = unname(group.mapping[f]) 
    meta.data[nrow(meta.data)+1, ] = c(f, this.group, unname(colors[this.group]))
    names(raw.data)[length(raw.data)] = f
  }
}


# Calculate Abundance and Diversity using package *alakazam*

## Define function to calculate abundance and diversity

is_outlier = function (x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

create_layout = function (start_num, end_num, num_col) {
  len_num = end_num - start_num + 1
  num_row = ceiling(len_num / num_col)
  tmp_vec = rep(NA, num_row * num_col)
  tmp_vec[1:len_num] = seq(start_num, end_num)
  return(matrix(tmp_vec, nrow = num_row, byrow = T))
}


performCalc = function (raw.data, meta.data, my.ylim=0.1, ds.num=200, extra.ylim=NA) {
  
  # plot cell count and cell type distribution
  cellnum.df = data.frame(sample = meta.data[, 'sample'], cellnum = unlist(lapply(raw.data, nrow)), stringsAsFactors=F)
  print(sprintf('Total %s cells', format(sum(cellnum.df$cellnum), big.mark=",", scientific=F)))
  
  # cell with none cell type will have a empty string as its cell type
  celltype.names = as.character(sort(unique(unlist(lapply(raw.data, function (x) unique(x[['annotation']]))))))
  
  my.data<-numeric()
  for (i in 1:length(raw.data)) {
    if(nrow(raw.data[[i]]) == 0) {
      temp.vect = rep(0, length(celltype.names))
      names(temp.vect) = celltype.names
    } else {
      temp.vect = table(as.character(raw.data[[i]]$annotation))
      temp.vect = temp.vect[celltype.names]
      temp.vect[is.na(temp.vect)] = 0
      names(temp.vect) = celltype.names
      temp.vect= t(as.matrix(temp.vect))[1, ]
    }
    temp.matrix = data.frame(cellnum = temp.vect, celltype = names(temp.vect), sample = rep(meta.data[i, 'sample'], length(temp.vect)))
    my.data = rbind(my.data, temp.matrix)
  }
  
  p1 = ggplot(data = my.data, aes(x = sample, y = cellnum, fill = celltype)) +
  geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(sample, cellnum, label = cellnum, fill = NULL), data = cellnum.df, vjust = -0.25, size = 3) + guides(fill = guide_legend(ncol = 1))
  
  # print('Plot of cell numbers')
  g = grid.arrange(p1, ncol=1)
  
  if (ds.num == 0) {return()}
  
  # calculate abundance and diversity
  abundance.list = list()
  diversity.list = list()
  
  for (i in 1:length(raw.data)) {
    
    if (nrow(raw.data[[i]]) < ds.num) {next}
    
    tmp = estimateAbundance(raw.data[[i]], ci=0.95, nboot=200, clone="clone", uniform=T, progress = F, group = NULL, min_n = ds.num, max_n = ds.num)
    diversity.list[[length(diversity.list)+1]] = alphaDiversity(raw.data[[i]], min_q=0, max_q=2, step_q=1, ci=0.95, nboot=200, clone="clone", uniform=T, group = NULL, min_n = ds.num, max_n = ds.num)
    
    # post-process abundance
    # exclude fake clones, re-ranking
    tmp.abundace = tmp@abundance %>%
      filter(!startsWith(clone, 'U'))
    tmp.abundace$rank = 1:nrow(tmp.abundace)
    
    tmp@abundance = tmp.abundace
    abundance.list[[length(abundance.list)+1]] = tmp
    
    names(abundance.list)[length(abundance.list)] = meta.data[i, 'sample']
    names(diversity.list)[length(diversity.list)] = meta.data[i, 'sample']
  }

  
  # plot abundance, each condition is a row
  p.list = list()
  if (!is.na(extra.ylim)) {p2.list = list()}
  
  num.col = 7
  layout.matrix = numeric()
  tmp = c(0)
  
  for (condition in names(colors)) {
    sample.names = meta.data[meta.data['group']==condition, 'sample']
    count = tmp[length(tmp)]
    for (sample.name in sample.names) {
      this.abundance = abundance.list[[sample.name]]
      if (!is.null(this.abundance)) {
        
        # get top 200 clones for ploting
        this.abundance@abundance = this.abundance@abundance[1:200, ]
        
        count = count + 1
        p.list[[length(p.list)+1]] = plotAbundanceCurve(this.abundance, colors=colors[condition], silent = TRUE) + ggtitle(sample.name) + ylim(0, my.ylim)
        if (!is.na(extra.ylim)) {
          p2.list[[length(p2.list)+1]] = plotAbundanceCurve(this.abundance, colors=colors[condition], silent = TRUE) + ggtitle(sample.name) + ylim(0, extra.ylim)
        }
      }
    }
    tmp = c(tmp, count)
    layout.matrix = rbind(layout.matrix, create_layout(tmp[length(tmp)-1]+1, tmp[length(tmp)], num.col))
  }

  # print('Plot abundance curves')
  g = grid.arrange(grobs = p.list, ncol=num.col, layout_matrix=layout.matrix)
  if (!is.na(extra.ylim)) {
    g2 = grid.arrange(grobs = p2.list, ncol=num.col, layout_matrix=layout.matrix)
  }
  
  
  # extract diversity
  diversity.matrix = numeric()
  for (i in 1:length(diversity.list)) {
    tmp = as.data.frame(diversity.list[[i]]@diversity)
    tmp[, 'group'] = unname(group.mapping[names(diversity.list)[i]])
    tmp[, 'sample'] = names(diversity.list)[i]
    diversity.matrix = rbind(diversity.matrix, tmp)
  }
  
  
  # plot diversity boxplot
  # and add outlier label
  index.names = c("Richness", "Shannon", "Simpson", "Shannon/Richness", "Simpson/Richness")
  p.list<-list()
  for (i in 1:3) {
    p.list[[i]] = diversity.matrix[diversity.matrix$q==i-1, ] %>%
      group_by(group) %>%
      mutate(outlier = ifelse(is_outlier(d), sample, NA)) %>%
      ggplot(., aes(x = group, y = d, color = group, label = sample)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(seed = 1)) +
      labs(x = 'Group', y = index.names[i]) +
      scale_color_manual(values = unname(colors)) +
      scale_x_discrete(limits = names(colors)) +
      theme(legend.position = "none") + 
      geom_text(position = position_jitter(seed = 1), hjust = -0.1, size = 3)
  }
  for (i in 4:5) {
    temp.matrix = diversity.matrix[diversity.matrix$q==i-3, ]
    temp.matrix.1 = diversity.matrix[diversity.matrix$q==0, ]
    temp.matrix[, 'd'] = temp.matrix[, 'd'] / temp.matrix.1[, 'd']
    p.list[[i]] = temp.matrix %>%
      group_by(group) %>%
      mutate(outlier = ifelse(is_outlier(d), sample, NA)) %>%
      ggplot(., aes(x = group, y = d, color = group, label = sample)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(seed = 1)) +
      labs(x = 'Group', y = index.names[i]) +
      scale_color_manual(values = unname(colors)) +
      scale_x_discrete(limits = names(colors)) +
      theme(legend.position = "none") +
      geom_text(position = position_jitter(seed = 1), hjust = -0.1, size = 3)
  }
  
  p.list[[6]] = cowplot::get_legend(ggplot(data=temp.matrix, aes(x = group, y = d, color = group, label = sample)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(seed = 1)) +
      labs(x = 'Group', y = index.names[i]) +
      scale_color_manual(values = unname(colors)) +
      scale_x_discrete(limits = names(colors)))
      
  #print('Boxplot of clone alpha diversity')
  do.call('grid.arrange', c(p.list, ncol = 3))
  
  return(list(abundance.list, diversity.matrix))
}


# NK and T cells from *subcluster_annotation3* (broad)
## Extract Ki67+ cells


nkt = read.table(file = file.path(home.dir, 'NKT_subcluster_cell_barcodes_091120_v3.csv'), header = TRUE, sep = ',', stringsAsFactors = F, row.names = 1)
nkt$cell_barcode = rownames(nkt)

sub.raw.data = list()
for (i in 1:nrow(meta.data)) {
  this.sample = meta.data[i, 'sample']
  this.nkt = nkt[nkt['orig.ident_1']==this.sample, ]
  if (meta.data[i, 'group'] == 'COVID19') {
    tmp.raw.data = raw.data[[i]]
    tmp.raw.data$match = tmp.raw.data$cell_id
    for (j in 1:nrow(tmp.raw.data)) {
      tmp.string = strsplit(tmp.raw.data[j, 'match'], '_')
      if (tmp.string[[1]][2] == 'nonCITE') {
        if (this.sample == 'TP8B') {
          tmp.raw.data[j, 'match'] = tmp.string[[1]][1]
        } else {
          tmp.raw.data[j, 'match'] = paste0(tmp.string[[1]][3], '_', tmp.string[[1]][1])
        }
      } else {
        tmp.sub.str = tmp.string[[1]][2]
        tmp.raw.data[j, 'match'] = paste0(substr(tmp.sub.str, 1, 4), '_', substr(tmp.sub.str, 5, 5), '_', tmp.string[[1]][1])
      }
    }
    tmp = merge(tmp.raw.data, this.nkt[, c('orig.barcode', 'cell_barcode', 'subcluster_annotation3')], by.x = 'match', by.y = 'orig.barcode')
  } else {
    this.nkt$match = paste0(this.nkt$orig.barcode, '_', this.sample)
    tmp = merge(raw.data[[i]], this.nkt[, c('match', 'cell_barcode', 'subcluster_annotation3')], by.x = 'cell_id', by.y = 'match')
  }
  tmp$annotation = tmp$subcluster_annotation3
  sub.raw.data[[length(sub.raw.data)+1]] = tmp
  names(sub.raw.data)[length(sub.raw.data)] = this.sample
}



## Prepare for output


getCellID = function (abundance, this.data, this.sample, this.condition, num=NA) {
  
  # get top num abundance clones, and add sample and condition infos
  if (is.na(num)) {
    this.abundance = abundance %>%
    mutate(sample = this.sample) %>%
    mutate(condition = this.condition) %>%
    subset(select = -c(group))
  } else {
    this.abundance = abundance[1:num, ] %>%
    mutate(sample = this.sample) %>%
    mutate(condition = this.condition) %>%
    subset(select = -c(group))
  }
  
  
  # outter join
  tmp = merge(x = this.abundance, y = this.data[, c('clone', 'cell_barcode')], by = 'clone', all.x = T)
  
  return(tmp[order(tmp$rank), ])
}

getAllCellID = function (abundance.list, data.list, this.condition, num=NA) {
  
  tmp.matrix = numeric()
  
  for (name in names(abundance.list)) {
    tmp = getCellID(abundance.list[[name]]@abundance, data.list[[name]], name, this.condition, num)
    tmp.matrix = rbind(tmp.matrix, tmp)
  }
  
  return(tmp.matrix)
}

# infos for output
cellid.matrix = numeric()
diversity.matrix = numeric()


## ki67 CD4 T cells

sub.cells = c('Ki67_CD4_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation==sub.cells), ])
. = performCalc(sub.data, meta.data, my.ylim = 0.25, ds.num = 0)


## ki67 CD8 T cells

sub.cells = c('Ki67_CD8_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation==sub.cells), ])
. = performCalc(sub.data, meta.data, my.ylim = 0.25, ds.num = 0)


## CD4+ T cells (excluding naive)

sub.cells = c('Memory_CD4_Tcells', 'Ki67_CD4_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation %in% sub.cells), ])
x = performCalc(sub.data, meta.data, my.ylim = 0.15, ds.num = 179, extra.ylim = 0.02)

# output diversity
this.condition = 'ki67+memory_CD4'
x[[2]]$condition = this.condition
diversity.matrix = rbind(diversity.matrix, x[[2]])

# output abundance
cellid.matrix = rbind(cellid.matrix, getAllCellID(x[[1]], sub.data, this.condition))


## CD8+ T cells (excluding naive)

sub.cells = c('Memory_CD8_Tcells', 'Ki67_CD8_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation %in% sub.cells), ])
x = performCalc(sub.data, meta.data, my.ylim = 0.45, ds.num = 123)

# output diversity
this.condition = 'ki67+memory_CD8'
x[[2]]$condition = this.condition
diversity.matrix = rbind(diversity.matrix, x[[2]])

# output abundance
cellid.matrix = rbind(cellid.matrix, getAllCellID(x[[1]], sub.data, this.condition))


## Naive CD4+ T cells

sub.cells = c('Naive_CD4_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation == sub.cells), ])
x = performCalc(sub.data, meta.data, my.ylim = 0.15, ds.num = 152, extra.ylim = 0.02)

# output diversity
this.condition = 'naive_CD4'
x[[2]]$condition = this.condition
diversity.matrix = rbind(diversity.matrix, x[[2]])

# output abundance
cellid.matrix = rbind(cellid.matrix, getAllCellID(x[[1]], sub.data, this.condition))


## Naive CD8+ T cells

sub.cells = c('Naive_CD8_Tcells')
sub.data = lapply(sub.raw.data, function (x) x[(!is.na(x$annotation)) & (x$annotation == sub.cells), ])
x = performCalc(sub.data, meta.data, my.ylim = 0.15, ds.num = 52, extra.ylim = 0.055)

# output diversity
this.condition = 'naive_CD8'
x[[2]]$condition = this.condition
diversity.matrix = rbind(diversity.matrix, x[[2]])

# output abundance
cellid.matrix = rbind(cellid.matrix, getAllCellID(x[[1]], sub.data, this.condition))


## Save results


write.csv(diversity.matrix, file = file.path(home.dir, 'diversity.csv'), row.names = F)

write.csv(cellid.matrix, file = file.path(home.dir, 'abundance_barcode.csv'), row.names = F)
