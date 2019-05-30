library(data.table)
library(stringr)
library(R.utils)


path <- 'data/ipsa_output/'
prefix <- 'output_'

ParseIpsaOuput <- function(feature, tissue, path, prefix, pull_with = 'sum'){
  dt.path <- paste0(path, '/', prefix, tissue, '_', feature, '.tsv')
  dt <- fread(dt.path)
  if(feature == 'exon'){
    dt$chr <- str_split_fixed(dt$exon_chr_start_end_strand, ':', 2)[, 1]
    dt$coord <- str_split_fixed(dt$exon_chr_start_end_strand, ':', 2)[, 2]
    dt$start <- as.numeric(str_split_fixed(dt$coord, '_', 2)[, 1])
    dt$end <- str_split_fixed(dt$coord, '_', 2)[, 2]
    dt$end <- as.numeric(unlist(lapply(dt$end, function(x){substr(x, 1, nchar(x) - 1)})))
    dt$strand <- unlist(lapply(dt$exon_chr_start_end_strand, function(x){substr(x, nchar(x), nchar(x))}))
    dt <- dt[, c('chr', 'start', 'end', 'strand', paste0(pull_with, '_exc_tumor'), paste0(pull_with, '_inc_tumor'), 
                 paste0(pull_with, '_exc_normal'), paste0(pull_with, '_inc_normal')), with = F]
    write.table(dt, paste0(path, '/', 'ipsa_', tissue, '_', feature, '_', pull_with, '.tsv'), sep = '\t', row.names = F, quote = F)
  }  else if(feature == 'site'){
    dt <- dt[, c('chr', 'coord', 'strand', 'sum_ret_tumor', 'sum_inc_tumor', 'sum_inc_normal', 'sum_ret_normal'), with = F]
    write.table(dt, paste0(path, '/', 'ipsa_', tissue, '_', feature, '.tsv'), sep = '\t', row.names = F, quote = F)
  } else if(feature == 'intron'){
    dt <- dt[, c('chr', 'start', 'end', 'strand', 'ret_tumor', 'inc_tumor', 'inc_normal', 'ret_normal'), with = F]
    write.table(dt, paste0(path, '/', 'ipsa_', tissue, '_', feature, '.tsv'), sep = '\t', row.names = F, quote = F)
  }
}

for (tissue in c('kidney', 'liver')){
  for (feature in c('site', 'exon', 'intron')){
    print(paste0(tissue, '_', feature))
    ParseIpsaOuput(feature, tissue, path, prefix, 'median')
  }
}

for (tissue in c('kidney', 'liver')){
  for (feature in c('exon')){
    print(paste0(tissue, '_', feature))
    ParseIpsaOuput(feature, tissue, path, prefix, 'median')
  }
}

