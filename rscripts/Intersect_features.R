#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(RColorBrewer)
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
  if(what %in% c('introns', 'circRNA')){
    dt$intersection.type <- 'crossing'
    if(nrow(dt[dt$percentage.ph == 1, ]) != 0) dt[dt$percentage.ph == 1, ]$intersection.type <- 'inside'
    if(nrow(dt[dt$percentage.interval == 1, ]) != 0) dt[dt$percentage.interval == 1, ]$intersection.type <- 'outside'    
  } else if(what %in% c('exons', 'transcript_starts', 'transcript_ends', 'polyA', 'CAGE', 'cryptic_ss')){
    dt$intersection.type <- 'crossing'
    if(nrow(dt[dt$percentage.interval == 1, ]) != 0) dt[dt$percentage.interval == 1, ]$intersection.type <- 'loop out'    
    #dt[dt$V13 == 0, ]$intersection.type <- 'not loop out'
    dt <- dt[dt$intersection.type == 'loop out', ]
  }
  if(what == 'cryptic_ss'){
    dt$intersection.type <- 'in handle'
  }
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


PlotTypes <- function(counts.all){
  tmp <- merge(subset(counts.all, experiment.type == 'real'), 
               subset(counts.all, experiment.type != 'real'), 
               by = 'intersection.type')
  tmp$ratio <- tmp$counts.x / tmp$counts.y
  g <- ggplot(tmp, aes(y = ratio, x = intersection.type)) + 
    geom_boxplot() +
    #geom_point(data = counts.all[counts.all$experiment.type == 'real', ], aes(y = counts), col = 'red') +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') +
    geom_hline(yintercept = 1, col = 'darkred', linetype = "dashed")
  g
}



DoIntesection <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db = 'counts'){
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
    if(what == 'cryptic_ss') genes <- conins 
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
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type, data = dt.all, FUN = sum))
  write.table(dt.all, paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_db_', db, '_position.tsv'), 
              sep = '\t', quote = F, row.names = F)
  p <- PlotTypes(counts.all)
  pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_', as.character(N.shuffle), '_db_', db, '.pdf'))
  print(p)
  dev.off()
}



##########################################
# For energy and spread cutoffs
##########################################

MakeCategories2 <- function(dt, what = 'energy', categories, strict = F){
  dt.new <- copy(dt)
  dt.new[, paste0(what, '.bin') := categories[1]]
  for(cat in categories[2:length(categories)]){
    if(strict){
      dt.tmp <- dt[dt[[what]] < cat, ]
    } else{
      dt.tmp <- dt[dt[[what]] <= cat, ]  
    }
    dt.tmp[, paste0(what, '.bin') := cat]
    dt.new <- rbind(dt.new, dt.tmp)
  }
  #dt.new <- dt.new[!duplicated(dt.new), ]
  dt.new
}


PlotTypesCutoffs <- function(counts.all, what. = 'introns', bin = 'energy.bin', colors. = energy.colors, shuffle){
  tmp <- merge(subset(counts.all, experiment.type == 'real'), 
               subset(counts.all, experiment.type != 'real'), 
               by = c('intersection.type', bin))
  tmp$ratio <- tmp$counts.x / tmp$counts.y
  title. <-  if(shuffle == 'phs') 'Random shift' else 'Random gene' 
  g <- ggplot(tmp, aes(y = ratio, x = intersection.type)) + 
    geom_boxplot(aes(col = get(bin))) +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') +
    ylab('enrichment') + 
    geom_hline(yintercept = 1, col = 'darkred', linetype = "dashed") +
    ggtitle(title.)
  if(bin == 'energy.bin'){
    g <- g + 
      scale_color_manual(values = colors., name = '\u0394G, kcal/mol')
  } else if(bin == 'spread.bim'){
    g <- g + 
      scale_color_manual(values = colors., name = 'spread, nts')
  } else if(bin == 'E-value.bin'){
    g <- g + 
      scale_color_manual(values = colors., name = 'E-value')
  }
  g
}

ProcessSS <- function(dt){
  path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
  phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
  phs.dt <- fread(phs)
  phs.dt$id <- str_split_fixed(phs.dt$V4, '_', 2)[, 1]
  phs.dt$energy <- as.numeric(str_split_fixed(phs.dt$V4, '_', 2)[, 2])
  dt$id <- as.character(dt$id)
  dt <- merge(dt, phs.dt[, c('id', 'energy'), with = F], by = 'id', all.x = T)
  evals.dt <- fread(evals)
  evals.dt$id <- as.character(evals.dt$id)
  dt <- merge(dt, evals.dt, by = 'id', all.x = T)
  dt$expression <- 'Low'
  dt[dt$counts >= 9, ]$expression <- 'High'
  dt
}



DoIntesectioncutoffs <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db = 'counts'){
  dt <- fread(paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_db_', db, '_position.tsv'))
  dt$id <- as.numeric(str_split_fixed(dt$ph.id, '_', 2)[, 1])
  dt$energy <- as.numeric(str_split_fixed(dt$ph.id, '_', 2)[, 2])
  
  phs.dt <- fread(phs)
  phs.dt$spread <- phs.dt$V3 - phs.dt$V2
  dt <- merge(dt, phs.dt[, c('V4', 'spread'), with = F], by.x = 'ph.id', by.y = 'V4', all.x = T)
  
  evals.dt <- fread(evals)
  dt <- merge(dt, evals.dt, by = 'id', all.x = T)
  #dt <- dt[dt$`E-value` <= 0.05, ]

  energy.cutoffs <- c(-15, -20, -25, -30)
  energy.colors <- brewer.pal(n = 12, "Paired")[c(4, 7, 8, 6)]
  names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')
  column = "energy"
  dt.energy <- MakeCategories2(dt = dt, what = column, categories = energy.cutoffs, strict = F)
  dt.energy$energy.bin <- paste('<=', dt.energy$energy.bin, sep = '')
  dt.energy$energy.bin <- factor(dt.energy$energy.bin)
  dt2 <- dt.energy[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'energy.bin'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + energy.bin, data = dt2, FUN = sum))
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_', N.shuffle, '_db_', db, '_cutoff_energy.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'energy.bin', colors. = energy.colors, shuffle)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_db_', db, '_energy_cutoffs', '.pdf'))
  print(p)
  dev.off()
  
  eval.cutoffs <- c(2, 1, 0.05)
  eval.colors <- brewer.pal(n = 9, "Set1")[c(9, 5, 2)]
  names(eval.colors) <- paste('<', eval.cutoffs, sep = '')
  names(eval.colors)[1] <- 'all'
  column = "E-value"
  dt.eval <- MakeCategories2(dt = dt, what = column, categories = eval.cutoffs, strict = T)
  dt.eval$`E-value.bin` <- paste('<', dt.eval$`E-value.bin`, sep = '')
  dt.eval[dt.eval$`E-value.bin` == '<2']$`E-value.bin` <- 'all'
  dt.eval$`E-value.bin` <- factor(dt.eval$`E-value.bin`)
  dt2 <- dt.eval[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'E-value.bin'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + `E-value.bin`, data = dt2, FUN = sum))
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_', N.shuffle, '_db_', db, '_cutoff_evalue.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'E-value.bin', colors. = eval.colors, shuffle)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_db_', db, '_evalue_cutoffs', '.pdf'))
  print(p)
  dev.off()

  
  # spread.cutoffs <- c(10000, 1000, 100)
  # spread.colors <- c('magenta', 'green', 'blue')
  # names(spread.colors) <- paste('<=', spread.cutoffs, sep = '')
  # column = "spread"
  # dt.spread <- MakeCategories2(dt = dt, what = column, categories = spread.cutoffs)
  # dt.spread$spread.bin <- paste('<=', dt.spread$spread.bin, sep = '')
  # dt.spread$spread.bin <- factor(dt.spread$spread.bin)
  # dt2 <- dt.spread[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'spread.bin'), with = F]
  # counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + spread.bin, data = dt2, FUN = sum))
  # write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_10_cutoff_spread.tsv'), sep = '\t', row.names = F)
  # p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'spread.bin', colors. = spread.colors)
  # cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_spread_cutoffs', '.pdf'))
  # print(p)
  # dev.off()
}


###### RUN
path.to.ph <- '../'
phs <- paste0(path.to.ph, '/pairs.bed')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'

path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'
conins <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.bed'
evals <- paste0(path.to.ph, 'R_scape_estended_all_pretty.tsv')

# path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
# phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered_inner.bed6')
# N.shuffle <- 10
# genes <- '../conservative_features/not_intersected_coding_genes.bed'


# exons
intervals <- '../conservative_features/gencode.v19.annotation_exons_coding.bed'
intervals <- '../exons_counts.bed'
intervals <- '../exons/exons_K562.bed'
# introns
intervals <- '../intron_counts.bed' # with counts ICGC
intervals <- '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed' # simple
intervals <- '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_shorter10000.bed' # simple short 
intervals <- '../introns/introns_K562.bed' # counts K562

db = 'simple'
what = 'introns'
for(shuffle in c('phs', 'genes')){
  intervals <- if(shuffle == 'phs'){
    '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_shorter10000.bed' 
  } else{
    '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed'
  }
  print(what)
  print(shuffle)
  DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
  DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
}


db = 'counts'
what = 'introns'
intervals <- '../introns/introns_K562_shorter10000.bed'
for(shuffle in c('phs', 'genes')){
  print(what)
  print(shuffle)
  print(db)
  DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
  DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
}

db = 'counts'
what = 'exons'
intervals <- '../exons/exons_K562.bed'
for(shuffle in c('phs', 'genes')){
  print(what)
  print(shuffle)
  DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
  DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
}

db = 'simple'
what = 'exons'
intervals <- '../conservative_features/gencode.v19.annotation_exons_coding.bed'
for(shuffle in c('phs', 'genes')){
  print(what)
  print(shuffle)
  DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
  DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
}

db = 'simple'
for(what in c('transcript_starts', 'transcript_ends')){
  if(what == 'transcript_starts'){
    intervals <- '../conservative_features/gencode.v19.annotation_transcripts_starts.bed'
  } else if(what == 'transcript_ends'){
    intervals <- '../conservative_features/gencode.v19.annotation_transcripts_ends.bed'
  }
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
    DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
  }
}

db = 'counts_20'
for(what in c('polyA')){
  if(what == 'polyA'){
    #intervals <- '../polyA/kidney_all.bed'
    intervals <- '../polyA/GSE30198_human.pas.uniq_filtered.bed'
  }
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
    DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
  }
}

db = 'counts'
for(what in c('CAGE')){
  if(what == 'CAGE'){
    intervals <- '../CAGE/NHEK_polyA+_all.bed'
  }
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins)
    DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)   
  }
}


what = 'circRNA'
for(db in c('CIRCpedia_K562', 'TSCD_liver')){
  if(db == 'CIRCpedia_K562'){
    intervals <- '../circRNA/CIRCpedia_K562.bed'
  } else{
    intervals <- '../circRNA/TSCD_liver.bed'
  }
  for(shuffle in c('phs', 'genes')){
    print(what)
    print(shuffle)
    print(db)
    DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins, db)
    DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db)   
  }  
}







path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered_handles.bed6')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'
conins <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.bed'
evals <- paste0(path.to.ph, 'R_scape_estended_all_pretty.tsv')

db = 'counts'
what = 'cryptic_ss'
intervals <- '../cryptic_ss/K562_cryptic_ss.bed'
shuffle = 'phs'
print(what)
print(shuffle)
DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins)
DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)   

db = 'simple'
what = 'cryptic_ss'
intervals <- '../cryptic_ss/cryptic_sites_GTEx.bed'
shuffle = 'phs'
print(what)
print(shuffle)
DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins)
DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)   


