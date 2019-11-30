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
  dt <- read.delim(pipe(paste0('bedtools intersect -a ', intervals, ' -b ', phs, ' -s -wao ')), header = F) 
  dt$ph.length <- dt$V9 - dt$V8 + 1
  dt$interval.length <- dt$V3 - dt$V2 + 1
  dt$percentage.interval <- dt$V13 / dt$interval.length 
  dt$intersection.type <- 'crossing'
  if(nrow(dt[dt$percentage.interval >= 0.9, ]) != 0) dt[dt$percentage.interval >= 0.9, ]$intersection.type <- 'looped out'    
  dt[dt$V8 == -1, ]$intersection.type <- 'not looped out'  
  dt <- subset(dt, intersection.type == 'looped out' | intersection.type == 'not looped out')
  dt
}



path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'
conins <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.bed'
evals <- paste0(path.to.ph, 'R_scape_estended_all_pretty.tsv')
intervals <- '../exons/exons_psi_K562.bed'


phs <- SelectWorkingSample(phs, genes)
print('selected phs')
intervals <- SelectWorkingSample(intervals, genes)
print('selected intervals')
dt <- IntersectionByType(phs, intervals, what = what)
print('intersected')
dt <- dt[order(dt$intersection.type, decreasing = F), ]
dt <- dt[!duplicated(dt[, c('V4')]), ]
g <- ggplot(dt, aes(x = V5, fill = intersection.type)) + 
  geom_histogram(aes(y=..density..), position = 'dodge') +
  xlab('PSI') +
  theme_linedraw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  theme(legend.title = element_blank()) +
  theme(legend.position = 'top')
wilcox.test(subset(dt, intersection.type == 'looped out')$V5, subset(dt, intersection.type != 'looped out')$V5, paired = F)

pdf('PSI_looped_out.pdf')
print(g)
dev.off()

