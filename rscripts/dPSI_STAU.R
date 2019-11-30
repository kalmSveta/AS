#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(Cairo)
options(scipen = 999)

SelectWorkingSample <- function(dt, genes){
  new.dt.name <- gsub('.bed', '_sample.bed', dt)
  system(paste0('bedtools intersect -a ', dt, ' -b ', genes, ' -f 0.9 -u -s > ', new.dt.name))
  return(new.dt.name)
}

IntersectionByType <- function(phs, intervals, what = 'introns'){
  dt <- read.delim(pipe(paste0('bedtools intersect -a ', phs, ' -b ', intervals, ' -s -wao ')), header = F) 
  dt$ph.length <- dt$V3 - dt$V2
  dt$interval.length <- dt$V9 - dt$V8
  dt <- dt[dt$V13 != 0 & dt$V13 != '.', ]
  dt$V13 <- as.numeric(as.character(dt$V13))
  dt$percentage.ph <- dt$V13 / dt$ph.length
  dt$percentage.interval <- dt$V13 / dt$interval.length
  if(nrow(dt[is.na(dt$percentage.interval), ]) != 0) dt[is.na(dt$percentage.interval), ]$percentage.interval <- 0
  if(what == 'introns'){
    dt$intersection.type <- 'crossing'
    if(nrow(dt[dt$percentage.ph == 1, ]) != 0) dt[dt$percentage.ph == 1, ]$intersection.type <- 'inside'
    if(nrow(dt[dt$percentage.interval == 1, ]) != 0) dt[dt$percentage.interval == 1, ]$intersection.type <- 'outside'    
  } else if(what %in% c('exons', 'transcript_starts', 'transcript_ends', 'polyA', 'CAGE', 'eClip')){
    dt$intersection.type <- 'crossing'
    if(nrow(dt[dt$percentage.interval == 1, ]) != 0) dt[dt$percentage.interval == 1, ]$intersection.type <- 'loop out'    
    #dt[dt$V13 == 0, ]$intersection.type <- 'not loop out'
    dt <- dt[dt$intersection.type == 'loop out', ]
  }
  if(what == 'eClip'){
    dt$intersection.type <- 'in handle'
  }
  dt
}

MakePretty <- function(dt, type){
  dt <- dt[, c('V4', 'V10', 'V11', 'intersection.type')]
  dt$experiment.type <- type
  colnames(dt) <- c('ph.id', 'interval.id', 'counts', 'intersection.type', 'experiment.type')
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
  dt$relative.end <- dt$V9 - dt$V3
  
  dt <- dt[dt$relative.end >= 0 & dt$relative.start >= 0, ]
  
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
  genes.dt[, new.gene.id := sample(V10), by = gene.bin]
  phs.dt <- merge(phs.dt, genes.dt[, c('V10', 'new.gene.id'), with = F], by = 'V10')
  phs.dt <- merge(phs.dt, genes.dt[, c('V10', 'new.gene.length', 'new.gene.start', 
                                       'new.gene.end', 'new.gene.chr', 'new.gene.strand'), with = F], 
                  by.x = 'new.gene.id', by.y = 'V10')
  
  # scale phs starts and ends to fit new gene length
  phs.dt$new.phs.start <- 0
  phs.dt$new.phs.end <- 0
  phs.dt$phs.length <- phs.dt$V3 - phs.dt$V2
  phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$new.phs.end <- 
    phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$new.gene.end - 
    phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$relative.end
  phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$new.phs.start <- 
    phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$new.phs.end - 
    phs.dt[phs.dt$relative.start >= phs.dt$relative.end, ]$phs.length
  
  phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$new.phs.start <- 
    phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$new.gene.start + 
    phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$relative.start
  phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$new.phs.end <- 
    phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$new.phs.start + 
    phs.dt[phs.dt$relative.start < phs.dt$relative.end, ]$phs.length
  
  phs.dt <- subset(phs.dt, new.phs.start >= new.gene.start)
  phs.dt <- subset(phs.dt, new.phs.end <= new.gene.end)
  
  # make bed format
  phs.dt <- phs.dt[, c('new.gene.chr', 'new.phs.start', 'new.phs.end', 'V4', 'V5', 'new.gene.strand'), with = F]
  phs.dt <- phs.dt[order(phs.dt$V4), ]
  phs.dt
}


PlotdPSi <- function(dt.all){
  g <- ggplot(dt.all, aes(y = counts, x = experiment.type)) + 
    geom_boxplot() +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') + 
    geom_hline(yintercept = 0, col = 'red')
  g
}


ggplot(dt.all, aes(x = counts, fill = experiment.type)) + 
  geom_histogram(bins = 100) +
  theme_linedraw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  xlab('') +
  xlim(-0.2, 0.2) +
  geom_vline(xintercept = 0)

DoIntesection <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph){
  phs <- SelectWorkingSample(phs, genes)
  print('selected phs')
  intervals <- SelectWorkingSample(intervals, genes)
  print('selected intervals')
  dt <- IntersectionByType(phs, intervals, what = what)
  print('intersected')
  dt.real <- MakePretty(dt, 'real')
  print('made pretty')
  
  dt.shuffled.list <- lapply(c(1:N.shuffle), function(iteration){
    print(iteration)
    shuffled.dt <- if(shuffle == 'phs') ShufflePhs(phs, genes, iteration) else ShuffleGenes(phs, genes, iteration)
    write.table(shuffled.dt, 'tmp.bed', sep = '\t', row.names = F, col.names = F, quote = F)
    dt <- IntersectionByType('tmp.bed', intervals, what = what)
    dt
  })
  
  dt.shuffled.list.pretty <- lapply(c(1:length(dt.shuffled.list)), function(i){
    dt <- dt.shuffled.list[[i]]
    dt.shuffled <- MakePretty(dt, as.character(i))
    dt.shuffled
  })
  list.all <- copy(dt.shuffled.list.pretty) 
  list.all[[length(list.all) + 1]] <- dt.real
  dt.all <- bind_rows(list.all)
  
  dt.all$energy <- as.numeric(str_split_fixed(dt.all$ph.id, '_', 2)[, 2])
  #dt.all <- dt.all[order(dt.all$energy), ]
  dt.all <- subset(dt.all, counts != 0)
  dt.all <- dt.all[!duplicated(dt.all[, c('interval.id', 'experiment.type')]), ]
  
  # dt.all$id <- str_split_fixed(dt.all$ph.id, '_', 2)[, 1]
  # qvalue <- fread('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/R_scape_estended_all_pretty.tsv')
  # dt.all <- merge(dt.all, qvalue, by = 'id', all.x = T)
  # dt.all <- dt.all[dt.all$`E-value` <= 0.01, ]
  # 
  # dt.all <- subset(dt.all, energy <= -25)
  
  write.table(dt.all, paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_position.tsv'), 
              sep = '\t', quote = F, row.names = F)
  
  p <- PlotdPSi(dt.all)
  pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '.pdf'))
  print(p)
  dev.off()
}


path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
N.shuffle <- 1
genes <- '../conservative_features/not_intersected_coding_genes.bed'


for(what in c('exons')){
  for(shuffle in c('phs', 'genes')){
    for(tissue in c('K562', 'HepG2')){
      intervals <- if(tissue == 'K562') '../STAU1/delta_psi_STAU1_K562.bed' else '../STAU1/delta_psi_STAU1_HepG2.bed'
      print(what)
      print(shuffle)
      print(tissue)
      DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)
    }
  }
}