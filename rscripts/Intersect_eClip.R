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
  dt$intersection.type <- 'crossing'
  if(nrow(dt[dt$percentage.interval >= 0.6, ]) != 0) dt[dt$percentage.interval >= 0.6, ]$intersection.type <- 'in handle'    
  dt <- dt[dt$intersection.type == 'in handle', ]
  dt
}

MakePretty <- function(dt, type){
  dt <- dt[, c('V4', 'V10', 'V11', 'intersection.type')]
  dt$experiment.type <- type
  colnames(dt) <- c('ph.id', 'intervals.id', 'counts', 'intersection.type', 'experiment.type')
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


PlotTypeseClip <- function(counts.all){
  tmp <- merge(subset(counts.all, experiment.type == 'real'), 
               subset(counts.all, experiment.type != 'real'), 
               by = c('intersection.type', 'gene'))
  tmp$ratio <- tmp$counts.x / tmp$counts.y
  medians <- aggregate(ratio ~ gene, data = tmp, FUN = median)
  medians <- medians[order(medians$ratio, decreasing = T), ]
  tmp$group <- 'overrepresented'
  tmp[tmp$gene %in% medians[medians$ratio < 1, ]$gene, ]$group <- 'underrepresented'
  tmp$group <- factor(tmp$group, levels = c('underrepresented', 'overrepresented'))
  tmp$gene <- factor(tmp$gene, levels = medians$gene)
  # to get pretty plots
  tmp1 <- tmp[tmp$gene %in% medians[medians$ratio < 1, ]$gene & tmp$ratio <= 1.2, ]
  tmp2 <- tmp[tmp$gene %in% medians[medians$ratio > 1, ]$gene & tmp$ratio <= 4.2, ]
  tmp <- rbind(tmp1, tmp2)
  #
  g <- ggplot(subset(tmp, p.values.x < 0.05), aes(y = ratio, x = gene)) + 
    facet_wrap(~group, nrow = 1, scales = 'free') + 
    geom_boxplot() +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') +
    ylab('enrichment') +
    geom_hline(yintercept = 1, col = 'darkred', linetype = "dashed") +
    theme(legend.position="top") +
    coord_flip()
  g
}

CalculatePvalues <- function(counts.all){
  p.values <- unlist(lapply(unique(counts.all$gene), function(gene.) 
  {
    if(length(unique(subset(counts.all, gene == gene. & experiment.type != 'real')$counts)) == 1 | 
       length(subset(counts.all, gene == gene. & experiment.type != 'real')$counts) == 0 |
       length(subset(counts.all, gene == gene. & experiment.type == 'real')$counts) == 0){
      return(NA)
    } else{
      t.test(subset(counts.all, gene == gene. & experiment.type != 'real')$counts, 
             mu = subset(counts.all, gene == gene. & experiment.type == 'real')$counts)$p.value  
    }
  }))
  names(p.values) <- unique(counts.all$gene)
  p.values <- p.adjust(p.values, method = 'BH')
  p.values <- as.data.frame(p.values)
  p.values$gene <- row.names(p.values) 
  counts.all <- merge(counts.all, p.values, by = 'gene')
  counts.all
}

# # add evals
# evals.dt <- fread(evals)
# dt.all$id <- str_split_fixed(dt.all$ph.id, '_', 2)[, 1]
# evals.dt$id <- as.character(evals.dt$id)
# dt.all <- merge(dt.all, evals.dt, by = 'id', all.x = T)
# dt.all <- dt.all[dt.all$`E-value` <= 0.05, ]

DoIntesection <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins){
  phs <- SelectWorkingSample(phs, genes)
  print('selected phs')
  intervals <- SelectWorkingSample(intervals, genes)
  print('selected intervals')
  conins <- SelectWorkingSample(conins, genes)
  print('selected conins')
  dt <- IntersectionByType(phs, intervals, what = what)
  print('intersected')
  dt.real <- MakePretty(dt, 'real')
  print('made pretty')
  
  dt.shuffled.list <- lapply(c(1:N.shuffle), function(iteration){
    print(iteration)
    shuffled.dt <- ShufflePhs(phs, conins, iteration)
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

  dt.all$gene <- str_split_fixed(dt.all$intervals.id, '_', 2)[, 1]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + gene, data = dt.all, FUN = sum))
  write.table(dt.all, paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_position.tsv'), 
              sep = '\t', quote = F, row.names = F)
  
  counts.all <- CalculatePvalues(counts.all)
  p <- PlotTypeseClip(counts.all)
  pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '.pdf'), width = 10, height = 10)
  print(p)
  dev.off()
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
  #dt.new <- dt.new[!duplicated(dt.new), ]
  dt.new
}

PlotTypesCutoffseCLIP <- function(counts.all, what. = 'introns', bin = 'energy.bin', colors. = energy.colors){
  tmp <- merge(subset(counts.all, experiment.type == 'real'), 
               subset(counts.all, experiment.type != 'real'), 
               by = c('intersection.type', bin, 'gene'))
  tmp$ratio <- tmp$counts.x / tmp$counts.y
  medians <- aggregate(ratio ~ gene, data = tmp, FUN = median)
  medians <- medians[order(medians$ratio, decreasing = T), ]
  tmp$gene <- factor(tmp$gene, levels = medians$gene)
  g <- ggplot(tmp, aes(y = ratio, x = gene, col = p.values.x < 0.05)) + 
    geom_boxplot() +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') +
    geom_hline(yintercept = 1, col = 'darkred', linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~energy.bin, nrow = 4)
  g
}


DoIntesectioncutoffseCLip <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph){
  dt <- fread(paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_position.tsv'))
  dt$id <- str_split_fixed(dt$ph.id, '_', 2)[, 1]
  phs.dt <- fread(phs)
  phs.dt$spread <- phs.dt$V3 - phs.dt$V2
  phs.dt$energy <- as.numeric(str_split_fixed(phs.dt$V4, '_', 2)[, 2])
  phs.dt$id <- str_split_fixed(phs.dt$V4, '_', 2)[, 1]
  dt <- merge(dt, phs.dt[, c('id', 'spread', 'energy'), with = F], by = 'id', all.x = T)
  
  energy.cutoffs <- c(-15, -20, -25, -30)
  energy.colors <- c('green', 'yellow', 'orange', 'red')
  names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')
  column = "energy"
  dt.energy <- MakeCategories2(dt = dt, what = column, categories = energy.cutoffs)
  dt.energy$energy.bin <- paste('<=', dt.energy$energy.bin, sep = '')
  dt.energy$energy.bin <- factor(dt.energy$energy.bin)
  dt2 <- dt.energy[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'energy.bin', 'gene'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + gene + energy.bin, 
                                        data = dt2, FUN = sum))
  tmp.list.w.pval <- lapply(unique(counts.all$energy.bin), function(bin) CalculatePvalues(subset(counts.all, energy.bin == bin)))
  counts.all <- Reduce(rbind, tmp.list.w.pval)
  
  
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_10_cutoff_energy.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffseCLIP(counts.all, what. = what, bin = 'energy.bin', colors. = energy.colors)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_energy_cutoffs', '.pdf'), width = 20, height = 12)
  print(p)
  dev.off()
  
  spread.cutoffs <- c(10000, 1000, 100)
  spread.colors <- c('magenta', 'green', 'blue')
  names(spread.colors) <- paste('<=', spread.cutoffs, sep = '')
  column = "spread"
  dt.spread <- MakeCategories2(dt = dt, what = column, categories = spread.cutoffs)
  dt.spread$spread.bin <- paste('<=', dt.spread$spread.bin, sep = '')
  dt.spread$spread.bin <- factor(dt.spread$spread.bin)
  dt2 <- dt.spread[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'spread.bin'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + spread.bin, data = dt2, FUN = sum))
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_10_cutoff_spread.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'spread.bin', colors. = spread.colors)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_spread_cutoffs', '.pdf'))
  print(p)
  dev.off()
}


path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered_handles.bed6')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'
conins <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.bed'
evals <- paste0(path.to.ph, 'R_scape_estended_all_pretty.tsv')

for(what in c('eClip')){
  if(what == 'eClip'){
    intervals <- '../eClip/peaks_merged_pretty.bed'
  }
  for(shuffle in c('phs')){
    print(what)
    print(shuffle)
    DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins)
    #DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)   
  }
}