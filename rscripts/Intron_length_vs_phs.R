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

intervals <- '../introns/introns_K562.bed' # counts K562

path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
genes <- '../conservative_features/not_intersected_coding_genes.bed'

phs <- SelectWorkingSample(phs, genes)
print('selected phs')
intervals <- SelectWorkingSample(intervals, genes)
print('selected intervals')
dt <- read.delim(pipe(paste0('bedtools intersect -a ', intervals, ' -b ', phs, ' -s -F 1 -wao ')), header = F) 
dt <- dt[!duplicated(dt[, c(1:6)]), ]
dt$intron.length <- dt$V3 - dt$V2
dt$ph <- F
dt[dt$V13 !=0, ]$ph <- T
dt <- dt[, c('V4', 'intron.length', 'ph')]
dt <- dt[order(dt$intron.length), ]
dt <- subset(dt, intron.length >= 20)
ggplot(dt, aes(x = log10(intron.length), fill = ph)) + 
  geom_histogram(aes(y=..density..), position = 'identity', alpha = 0.5,  bins = 50)

dt.ph <- subset(dt, ph)
dt.noph <- subset(dt, ph == F)
