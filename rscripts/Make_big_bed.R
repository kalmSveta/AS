#!/usr/bin/Rscript
library(data.table)
library(stringr)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)

# Usage: Rscript Make_big_bed.R path.to.ph.main track_name name_columns, path.to.chr

path.to.ph.main <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.tsv'
path.to.chr <- '../tools/hg19.chrom.sizes'
track_name <- 'filtered panhandles'
name_columns <- 'id,energy'

path.to.ph.main <- args[1]
track_name <- args[2]
name_columns <- args[3]
path.to.chr <- args[4]

name_columns <- unlist(strsplit(name_columns, ',', fixed = T))

Make_bed = function(dt, file_name, track_name, name_columns){
  dt[, index:=do.call(paste, c(.SD, sep = "_")), .SDcols = name_columns]
  dt$score <- '0'
  dt$block_count <- 2
  dt$blockSizes <- paste(as.character(dt$al1_length), as.character(dt$al2_length - 1), sep = ',')
  dt$blockStarts <- paste(as.character(0),as.character(dt$panhandle_right_hand - dt$panhandle_start), sep = ',')
  if(! 'itemRGB' %in% colnames(dt)){
    dt$itemRGB <- '0,0,0'
  }
  bed <- dt[, c('chr', 'panhandle_start', 'panhandle_end', 'index', 'score', 'strand', 'panhandle_start', 'panhandle_end', 'itemRGB',
                'block_count', 'blockSizes', 'blockStarts')]
  f <- file(file_name, "w")
  writeLines(paste("track name= ", track_name, sep = ''), f)
  write.table(bed, f, sep = '\t', row.names = F, col.names = F, quote = F)  
  close(f)
  return(0)
}

BedToBigBed <- function(path.to.ph.bed, n_header_lines = 0, path.to.chr){
  x <- paste0('tail -n +', n_header_lines + 1, ' ', path.to.ph.bed, 
              ' | bedtools sort -i stdin > sorted.bed')
  system(x)
  x <- paste0('../tools/bedToBigBed sorted.bed ', path.to.chr, ' ', gsub('bed', 'bb', path.to.ph.bed))
  system(x)
  system('rm sorted.bed')
}

dt <- fread(path.to.ph.main)
Make_bed(dt, gsub('tsv', 'bed', path.to.ph.main), track_name, name_columns)
BedToBigBed(gsub('tsv', 'bed', path.to.ph.main), n_header_lines = 1, path.to.chr)






