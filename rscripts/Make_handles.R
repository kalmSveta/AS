#!/usr/bin/Rscript
library(data.table)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)

# Make_handles.R path.to.ph.main 
path.to.ph.main <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.tsv'
path.to.ph.main <- args[1]

MakeHandlesBed <- function(dt){
  dt.left <- dt[, c('chr', 'panhandle_start', 'panhandle_left_hand', 'id')]
  dt.right <- dt[, c('chr', 'panhandle_right_hand', 'panhandle_end', 'id')]
  colnames(dt.left) <- c('chr', 'start', 'stop', 'name')
  colnames(dt.right) <- c('chr', 'start', 'stop', 'name')
  dt.left$name <- paste0(dt.left$name, '_left')
  dt.right$name <- paste0(dt.right$name, '_right')
  dt.both <- rbind(dt.left, dt.right)
  dt.both
}

dt <- fread(path.to.ph.main)
handles <- MakeHandlesBed(dt)
write.table(handles, gsub('.tsv', '_handles.bed', path.to.ph.main), sep = '\t', row.names = F, col.names = F, quote = F)
