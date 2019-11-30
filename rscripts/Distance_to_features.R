#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(plyr)
library(Cairo)
library(Hmisc)
options(scipen = 999)

# args <- commandArgs(trailingOnly = TRUE)
# 
# phs <-  args[1]
# intervals <- args[2]
# genes <- args[3]
# N.shuffle <- args[4]
# what <- args[5]
# shuffle <- args[6]


# select phs and intervals in genes
SelectWorkingSample <- function(dt, genes){
  new.dt.name <- gsub('.bed', '_sample.bed', dt)
  system(paste0('bedtools intersect -a ', dt, ' -b ', genes, ' -f 1 -u -s > ', new.dt.name))
  return(new.dt.name)
}

CountDistancesIntrons <- function(phs, intervals){
  dt <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', phs, ' -b ', intervals, ' -s -wo -f 1')), header = F))
  left <- aggregate(V8 ~ V4, data = dt, FUN = max)
  colnames(left) <- c('ph.id', 'left.coord')
  right <- aggregate(V9 ~ V4, data = dt, FUN = min)
  colnames(right) <- c('ph.id', 'right.coord')
  coords <- merge(left, right, by = 'ph.id')
  dt <- merge(coords, dt, by.x = 'ph.id', by.y = 'V4')
  dt$left.distance <- dt$V2 - dt$left.coord
  dt$right.distance <- dt$right.coord - dt$V3
  dt <- dt[dt$left.coord == dt$V8 & dt$right.coord == dt$V9, ]
  # correct on strand
  dt$five.prime.distance <- dt$left.distance
  dt[dt$V6  == '-', ]$five.prime.distance <- dt[dt$V6  == '-', ]$right.distance
  dt$three.prime.distance <- dt$right.distance
  dt[dt$V6  == '-', ]$three.prime.distance <- dt[dt$V6  == '-', ]$left.distance
  dt <- dt[, c('ph.id', 'five.prime.distance', 'three.prime.distance', 'V11')]
  colnames(dt)[4] <- 'counts' 
  dt
}

CountDistancesExons <- function(phs, intervals){
  dt <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', phs, ' -b ', intervals, ' -s -wo -F 1')), header = F))
  left <- aggregate(V8 ~ V4, data = dt, FUN = min)
  colnames(left) <- c('ph.id', 'left.coord')
  right <- aggregate(V9 ~ V4, data = dt, FUN = max)
  colnames(right) <- c('ph.id', 'right.coord')
  coords <- merge(left, right, by = 'ph.id')
  dt <- merge(coords, dt, by.x = 'ph.id', by.y = 'V4')
  dt$left.distance <- dt$left.coord - dt$V2 
  dt$right.distance <- dt$V3 - dt$right.coord
  dt <- dt[dt$left.coord == dt$V8 & dt$right.coord == dt$V9, ]
  # correct on strand
  dt$five.prime.distance <- dt$left.distance
  dt[dt$V6  == '-', ]$five.prime.distance <- dt[dt$V6  == '-', ]$right.distance
  dt$three.prime.distance <- dt$right.distance
  dt[dt$V6  == '-', ]$three.prime.distance <- dt[dt$V6  == '-', ]$left.distance
  dt <- dt[, c('ph.id', 'five.prime.distance', 'three.prime.distance', 'V11')]
  colnames(dt)[4] <- 'counts' 
  dt <- dt[! duplicated(dt), ]
  dt
}

CountDistanceseClip <- function(phs, intervals){
  system(paste0('bedtools sort -i ', phs, ' > tmp'))
  system(paste0('mv tmp ', phs))
  dt <- as.data.table(read.delim(pipe(paste0('bedtools sort -i ', intervals, '| bedtools closest -a ', phs, ' -b stdin', ' -s -io -d -t first')), header = F))
  dt <- dt[, c('V4', 'V13', 'V11'), with = F]
  colnames(dt) <- c('ph.id', 'distance', 'counts')
  dt
}

ShufflePhs <- function(phs, genes, seed){
  set.seed(seed)
  dt <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', phs, ' -b ', genes, ' -s -wa -wb ')), header = F))
  dt$ph.length <- dt$V3 - dt$V2
  dt$random <- runif(nrow(dt), -1, 1)
  dt$random.scaled <- dt$random
  dt[dt$random >= 0, ]$random.scaled <- dt[dt$random >= 0, ]$random.scaled * (dt[dt$random >= 0, ]$V9 - dt[dt$random >= 0, ]$V3)
  dt[dt$random < 0, ]$random.scaled <- dt[dt$random < 0, ]$random * (dt[dt$random < 0, ]$V2 - dt[dt$random < 0, ]$V8)
  dt$random.scaled <- round(dt$random.scaled)
  dt$new.start <- dt$V2 + dt$random.scaled
  dt$new.end <- dt$V3 + dt$random.scaled
  dt <- dt[, c('V1', 'new.start', 'new.end', 'V4', 'V5', 'V6'), with = F]
  dt
}

ShuffleGenes <- function(phs, genes, seed, n.genes.in.bin = 10){
  set.seed(seed)
  # fin gene for every ph
  dt <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', phs, ' -b ', genes, ' -s -wa -wb ')), header = F))
  dt$old.gene.length <- dt$V9 - dt$V8
  
  # divide genes in bins by their length
  x <- copy(dt)
  x <- x[!duplicated(dt$V10), ]
  x <- x[order(x$old.gene.length), ]
  n.bins <- round(length(x$old.gene.length) / n.genes.in.bin) 
  x$gene.bin <- as.numeric(cut_number(x$old.gene.length, n = n.bins, labels = 1:n.bins))
  dt <- merge(dt, x[, c('V10', 'gene.bin'), with = F], by = 'V10')
  
  # find old relative start and length of phs
  dt$relative.start <- dt$V2 - dt$V8
  dt$relative.end <- dt$V3 - dt$V8
  
  # make genes dt and phs dt
  genes.dt <- dt[, c('V7', 'V8', 'V9', 'V10', 'V11', 'V12', 'gene.bin', 'old.gene.length'), with = F]
  colnames(genes.dt)[colnames(genes.dt) == 'old.gene.length'] <- 'new.gene.length'
  colnames(genes.dt)[colnames(genes.dt) == 'V8'] <- 'new.gene.start'
  colnames(genes.dt)[colnames(genes.dt) == 'V9'] <- 'new.gene.end'
  colnames(genes.dt)[colnames(genes.dt) == 'V12'] <- 'new.gene.strand'
  colnames(genes.dt)[colnames(genes.dt) == 'V7'] <- 'new.gene.chr'
  genes.dt <- genes.dt[!duplicated(genes.dt$V10), ]
  phs.dt <- dt[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'old.gene.length', 'relative.start', 'relative.end', 'V10', 'gene.bin'), with = F]
  
  # shuffle gene ids in phs dt
  phs.dt[, new.gene.id := sample(V10), by = gene.bin]
  phs.dt <- merge(phs.dt, genes.dt[, c('V10', 'new.gene.length', 'new.gene.start', 
                                       'new.gene.end', 'new.gene.chr', 'new.gene.strand'), with = F], 
                  by.x = 'new.gene.id', by.y = 'V10')
  
  # scale phs starts and ends to fit new gene length
  phs.dt$new.phs.start <- round(phs.dt$new.gene.start + (phs.dt$new.gene.length / phs.dt$old.gene.length) * phs.dt$relative.start)
  phs.dt$new.phs.end <- round(phs.dt$new.gene.start + (phs.dt$new.gene.length / phs.dt$old.gene.length) * phs.dt$relative.end)
  
  # make bed format
  phs.dt <- phs.dt[, c('new.gene.chr', 'new.phs.start', 'new.phs.end', 'V4', 'V5', 'new.gene.strand'), with = F]
  phs.dt <- phs.dt[order(phs.dt$V4), ]
  phs.dt
}


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


PlotTypes <- function(dt.melt, title){
  dt.melt.five <- dt.melt[dt.melt$end == 'five.prime.distance', ]
  dt.melt.three <- dt.melt[dt.melt$end == 'three.prime.distance', ]
  
  by.type.five <- ddply(dt.melt.five, c("type"), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                          probs=c(0, .25, .5, .75, 1)))
  colnames(by.type.five) <- gsub(' ', '', colnames(by.type.five))
  by.type.three <- ddply(dt.melt.three, c("type"), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                            probs=c(0, .25, .5, .75, 1)))
  colnames(by.type.three) <- gsub(' ', '', colnames(by.type.three))
  
  
  p5 <- ggplot(by.type.five, aes(x = type, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                 upper = `75%`,  ymax = `100%`, fill = type)) + 
    geom_boxplot(stat = "identity") +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_fill_discrete(name = "", labels = c('expected', 'observed')) +
    #scale_fill_discrete(name = "") +
    ylab('distance, nts') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_flip()
  
  p3 <- ggplot(by.type.three, aes(x = type, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                  upper = `75%`,  ymax = `100%`, fill = type)) + 
    geom_boxplot(stat = "identity") + 
    scale_y_reverse(breaks = seq(0, 10000, 2500)) +
    scale_fill_discrete(name = "", labels = c('expected', 'observed')) +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab('distance, nts') +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.position = "none")
  
  legend <- get_legend(p5)
  p5 <- p5 + theme(legend.position="none")
  
  
  pdf(gsub('/', '.', paste0('distance_', what, '_shuffle_', shuffle, '_', title, '.pdf')))
  grid.arrange(p5, p3, legend, nrow = 1, widths=c(2.3, 2.3, 0.8), 
               bottom = textGrob(title, gp = gpar(fontsize = 20, font = 3)))
  dev.off()
}

PlotTypeseClip <- function(dt.melt, title){
  by.type <- ddply(dt.melt, c("type"), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                          probs=c(0, .25, .5, .75, 1)))
  colnames(by.type) <- gsub(' ', '', colnames(by.type))
  p <- ggplot(by.type, aes(x = type, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                 upper = `75%`,  ymax = `100%`, fill = type)) + 
    geom_boxplot(stat = "identity") +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_fill_discrete(name = "", labels = c('expected', 'observed')) +
    ylab('distance, nts') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_flip(ylim = c(0, 5000))
  
  pdf(gsub('/', '.', paste0('distance_', what, '_shuffle_', shuffle, '_', title, '.pdf')))
  print(p)
  dev.off()
}

DoDistance <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph){
  phs <- SelectWorkingSample(phs, genes)
  intervals <- SelectWorkingSample(intervals, genes)
  distances.real <- 
    if(what == 'introns') {
    CountDistancesIntrons(phs, intervals) 
  } else if(what == 'exons'){
    CountDistancesExons(phs, intervals)
  } else if(what == 'eClip'){
    CountDistanceseClip(phs, intervals)
  }
  distances.real$type <- 'real'
  
  dt.shuffled.list <- lapply(c(1:N.shuffle), function(iteration){
    print(iteration)
    shuffled.dt <- if(shuffle == 'phs') ShufflePhs(phs, genes, iteration) else ShuffleGenes(phs, genes, iteration)
    write.table(shuffled.dt, 'tmp.bed', sep = '\t', row.names = F, col.names = F, quote = F)
    dt <-     
      if(what == 'introns') {
      CountDistancesIntrons('tmp.bed', intervals) 
    } else if(what == 'exons'){
      CountDistancesExons('tmp.bed', intervals)
    } else if(what == 'eClip'){
      CountDistanceseClip('tmp.bed', intervals)
    }
    dt$type <- as.character(iteration)
    dt
  })
  print('shuffled')
  list.all <- copy(dt.shuffled.list)
  list.all[[length(list.all) + 1]] <- distances.real
  dt.all <- bind_rows(list.all)
  write.table(dt.all, paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_distance.tsv'),
             sep = '\t', quote = F, row.names = F)
  print('wrote.table')
  if(what == 'eClip'){
    dt.all.melt <- copy(dt.all)
    dt.all.melt <- dt.all.melt[dt.all.melt$distance <= 10000, ]
    PlotTypeseClip(dt.all.melt, 'general')
    print('plotted')
    w <- wilcox.test(dt.all.melt[dt.all.melt$type == '1', ]$distance, 
                     dt.all.melt[dt.all.melt$type == 'real', ]$distance, paired = F)
    print(w)
  } else{
    dt.all.melt <- melt(dt.all, id.vars = c('ph.id', 'counts', 'type'), 
                        measure.vars = c('five.prime.distance', 'three.prime.distance'), 
                        variable.name = 'end')    
    colnames(dt.all.melt)[colnames(dt.all.melt) == 'value'] <- 'distance'
    dt.all.melt <- dt.all.melt[dt.all.melt$distance <= 10000, ]
    PlotTypes(dt.all.melt, 'general')
    print('plotted')
    w <- wilcox.test(dt.all.melt[dt.all.melt$type == '1' & dt.all.melt$end == 'five.prime.distance', ]$distance, 
                     dt.all.melt[dt.all.melt$type == 'real' & dt.all.melt$end == 'five.prime.distance', ]$distance, paired = F)
    print(w)
    w <- wilcox.test(dt.all.melt[dt.all.melt$type == '1' & dt.all.melt$end == 'three.prime.distance', ]$distance, 
                     dt.all.melt[dt.all.melt$type == 'real' & dt.all.melt$end == 'three.prime.distance', ]$distance, paired = F)  
    print(w)
  }
}

##########################################
# For energy and spread cutoffs
##########################################
MakeCategories2 <- function(dt, what = 'energy', categories){
  dt.new <- copy(dt)
  dt.new[, paste0(what, '.bin') := categories[1]]
  for(cat in categories[2:length(categories)]){
    dt.tmp <- dt.new[dt.new[[what]] <= cat, ]
    dt.tmp[, paste0(what, '.bin') := cat]
    dt.new <- rbind(dt.new, dt.tmp)
  }
  dt.new
}


PlotTypesCutoffs <- function(dt.melt, title, column, colors.){
  dt.melt.five <- dt.melt[dt.melt$end == 'five.prime.distance', ]
  dt.melt.three <- dt.melt[dt.melt$end == 'three.prime.distance', ]

  
  by.type.five <- ddply(dt.melt.five, c("type", column), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                          probs=c(0, .25, .5, .75, 1)))
  colnames(by.type.five) <- gsub(' ', '', colnames(by.type.five))
  by.type.three <- ddply(dt.melt.three, c("type", column), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                            probs=c(0, .25, .5, .75, 1)))
  colnames(by.type.three) <- gsub(' ', '', colnames(by.type.three))
  
  by.type.five$type <- factor(by.type.five$type)
  by.type.five[, column] <- factor(by.type.five[, column])
  colnames(by.type.five)[colnames(by.type.five) == column] <- 'bin'
  by.type.three$type <- factor(by.type.three$type)
  by.type.three[, column] <- factor(by.type.three[, column])
  colnames(by.type.three)[colnames(by.type.three) == column] <- 'bin'
  
  p5 <- ggplot(by.type.five, aes(x = bin, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                 upper = `75%`,  ymax = `100%`, fill = bin, color = type)) + 
    geom_boxplot(stat = "identity") +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_color_discrete(name = "", labels = c('expected', 'observed')) +
    ylab('distance, nts') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_flip(ylim = c(0, 5500)) + 
    guides(fill = guide_legend(reverse = TRUE))
  
  if(column == 'energy.bin'){
    p5 <- p5 + 
      scale_fill_manual(name = '\u0394G, kcal/mol', values = colors.)
  } else{
    p5 <- p5 + 
      scale_fill_manual(name = 'spread, nts', values = colors.)
  }
  
  p3 <- ggplot(by.type.three, aes(x = bin, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                  upper = `75%`,  ymax = `100%`, fill = bin, color = type)) + 
    geom_boxplot(stat = "identity", lwd = 1.5) + 
    scale_y_reverse(breaks = seq(0, 10000, 2500)) +
    scale_fill_manual(name = "", values = colors.) +
    scale_color_discrete(name = "", labels = c('expected', 'observed')) +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab('distance, nts') +
    coord_flip(ylim = c(0, 5500)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.position = "none")
  
  legend <- get_legend(p5)
  p5 <- p5 + theme(legend.position="none") + geom_boxplot(stat = "identity", lwd = 1.5)
  
  
  cairo_pdf(gsub('/', '.', paste0('distance_', what, '_shuffle_', shuffle, '_', title, '.pdf')))
  grid.arrange(p5, p3, legend, nrow = 1, widths=c(2.3, 2.3, 0.8), 
               bottom = textGrob(title, gp = gpar(fontsize = 20, font = 3)))
  dev.off()
}

PlotTypesCutoffseClip <- function(dt.melt, title, column, colors.){
  by.type.five <- ddply(dt.melt, c("type", column), function(x) wtd.quantile(x$distance, weights=x$counts, 
                                                                                  probs=c(0, .25, .5, .75, 1)))
  colnames(by.type.five) <- gsub(' ', '', colnames(by.type.five))
  by.type.five$type <- factor(by.type.five$type)
  by.type.five[, column] <- factor(by.type.five[, column])
  colnames(by.type.five)[colnames(by.type.five) == column] <- 'bin'
  
  p <- ggplot(by.type.five, aes(x = bin, ymin = `0%`, lower = `25%`, middle = `50%`, 
                                 upper = `75%`,  ymax = `100%`, fill = bin, color = type)) + 
    geom_boxplot(stat = "identity") +
    theme_linedraw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_color_discrete(name = "", labels = c('expected', 'observed')) +
    ylab('distance, nts') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_flip(ylim = c(0, 3000)) + 
    guides(fill = guide_legend(reverse = TRUE))
  
  if(column == 'energy.bin'){
    p <- p + 
      scale_fill_manual(name = '\u0394G, kcal/mol', values = colors.)
  } else{
    p <- p + 
      scale_fill_manual(name = 'spread, nts', values = colors.)
  }
  cairo_pdf(gsub('/', '.', paste0('distance_', what, '_shuffle_', shuffle, '_', title, '.pdf')))
  print(p)
  dev.off()
}

DoDistancecutoffs <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph){
  dt <- fread(paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_distance.tsv'))
  dt$id <- as.numeric(str_split_fixed(dt$ph.id, '_', 2)[, 1])
  dt$energy <- as.numeric(str_split_fixed(dt$ph.id, '_', 2)[, 2])
  phs.dt <- fread(phs)
  phs.dt$spread <- phs.dt$V3 - phs.dt$V2
  dt <- merge(dt, phs.dt[, c('V4', 'spread'), with = F], by.x = 'ph.id', by.y = 'V4', .x = T)
  dt.melt <- melt(dt, id.vars = c('ph.id', 'counts', 'type', 'energy', 'spread'), 
                  measure.vars = c('five.prime.distance', 'three.prime.distance'), 
                  variable.name = 'end')
  colnames(dt.melt)[colnames(dt.melt) == 'value'] <- 'distance'
  dt.melt$end <- factor(dt.melt$end)
  write.table(dt.melt, paste0('distance_', what, '_shuffle_', shuffle, '_cutoffs.tsv'), sep = '\t', row.names = F)
  
  dt.melt <- dt.melt[dt.melt$distance <= 10000, ]
  
  energy.cutoffs <- c(-15, -20, -25, -30)
  energy.colors <- c('green', 'yellow', 'orange', 'red')
  names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')
  dt.melt.energy <- MakeCategories2(dt.melt, what = 'energy', energy.cutoffs)
  dt.melt.energy$energy.bin <- paste('<=', dt.melt.energy$energy.bin, sep = '')
  PlotTypesCutoffs(dt.melt.energy, 'spread <= 10.000 nts', 'energy.bin', energy.colors)
  
  spread.cutoffs <- c(10000, 1000, 100)
  spread.colors <- c('magenta', 'green', 'blue')
  names(spread.colors) <- paste('<=', spread.cutoffs, sep = '')
  dt.melt.spread <- MakeCategories2(dt.melt, what = 'spread', spread.cutoffs)
  dt.melt.spread$spread.bin <- paste('<=', dt.melt.spread$spread.bin, sep = '')
  PlotTypesCutoffs(dt.melt.spread, 'energy <= -15kcal/mol', 'spread.bin', spread.colors)
}


############ RUN
path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
N.shuffle <- 1
genes <- '../conservative_features/not_intersected_coding_genes.bed'

path.to.ph <- '../'
phs <- paste0(path.to.ph, '/pairs.bed')

for(what in c('introns')){
  if(what == 'introns'){
    intervals <- '../intron_counts.bed'
    #intervals <- '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed'
  } else if(what == 'exons'){
    intervals <- '../exon_counts.bed'
  }
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    DoDistance(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)
    DoDistancecutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)
  }
}
path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered_handles.bed6')
N.shuffle <- 1
genes <- '../conservative_features/not_intersected_coding_genes.bed'

for(what in c('eClip')){
  intervals <- '../eClip/peaks_merged_pretty.bed'
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    DoDistance(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)
    #DoDistancecutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)
  }
}

