#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(stringr)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)

# Usage: Make_big_bed.R path.to.ph.main name.columns color.energy energy.threshold path.to.chr

# path.to.ph.main <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.tsv'
# path.to.chr <- '../tools/hg19.chrom.sizes'
# name.columns <- 'id,energy'
# energy.threshold <- -15
# color.energy <- T

path.to.ph.main <- args[1]
name.columns <- args[2]
color.energy <- as.logical(args[3])
energy.threshold <- as.numeric(args[4])
path.to.chr <- args[5]


name.columns <- unlist(strsplit(name.columns, ',', fixed = T))

ColorPh <- function(dt){
  color.names <- c('green', 'yellow', 'orange', 'red', 'magenta')
  color.values <- c('0,100,0', '255,200,0', '255,100,0', '255,0,0', '255,0,255')
  names(color.values) <- color.names
  dt$itemRGB <- color.values['green']
  dt[dt$energy <= -20 & dt$energy > -25, ]$itemRGB <- color.values['yellow']
  dt[dt$energy <= -25 & dt$energy > -30, ]$itemRGB <- color.values['orange']
  dt[dt$energy <= -30 & dt$energy > -35, ]$itemRGB <- color.values['red']
  dt[dt$energy <= -35, ]$itemRGB <- color.values['magenta']
  dt
}

MakeBed = function(dt, file_name, name.columns, energy.threshold){
  dt[, index:=do.call(paste, c(.SD, sep = "_")), .SDcols = name.columns]
  dt$score <- '0'
  dt$block_count <- 2
  dt$blockSizes <- paste(as.character(dt$al1_length), as.character(dt$al2_length - 1), sep = ',')
  dt$blockStarts <- paste(as.character(0),as.character(dt$panhandle_right_hand - dt$panhandle_start), sep = ',')
  dt <- dt[dt$energy <= energy.threshold, ]
  if(! 'itemRGB' %in% colnames(dt)){
    dt$itemRGB <- '0,0,0'
  }
  bed <- dt[, c('chr', 'panhandle_start', 'panhandle_end', 'index', 'score', 'strand', 'panhandle_start', 'panhandle_end', 'itemRGB',
                'block_count', 'blockSizes', 'blockStarts')]
  f <- file(file_name, "w")
  #writeLines(paste("track name= ", track_name, sep = ''), f)
  write.table(bed, f, sep = '\t', row.names = F, col.names = F, quote = F)  
  close(f)
}

BedToBigBed <- function(path.to.ph.bed, n.header.lines = 0, path.to.chr){
  x <- paste0('tail -n +', n.header.lines + 1, ' ', path.to.ph.bed, 
              ' | bedtools sort -i stdin > sorted.bed')
  system(x)
  x <- paste0('../tools/bedToBigBed sorted.bed ', path.to.chr, ' ', gsub('bed', 'bb', path.to.ph.bed))
  system(x)
  system('rm sorted.bed')
}

dt <- fread(path.to.ph.main)
if(color.energy) ColorPh(dt) -> dt
out.file.name <- gsub('.tsv', paste0('_dGcutoff', energy.threshold, '.bed'), path.to.ph.main)
MakeBed(dt, out.file.name, name.columns, energy.threshold)
BedToBigBed(out.file.name, n.header.lines = 0, path.to.chr)






