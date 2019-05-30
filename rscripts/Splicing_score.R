library(data.table)
library(stringr)
options(scipen = 999)

path.to.ph <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed.tsv'
ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered_handles.bed'
ph.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered.bed'
ph.handle.flanked.path.bed <- gsub('.bed', '_with_flanks.bed', ph.handle.path.bed)
prefix <- 'data/ipsa_output/ipsa_'

exons.bridges.path <- '../RBPs/eCLIP_and_KD.tsv'
ph <- fread(path.to.ph)


LoopOut <- function(tissue, feature, prefix){
  feature.path <- paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')
  feature.dt <- fread(feature.path)
  x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
  system(x)  
  overlapped <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
  ncol <- dim(feature.dt)[2]
  overlapped$overlap <- T
  if(!any(grepl('ph.id', colnames(overlapped)))){
    overlapped$ph.id <- as.numeric(str_split_fixed(overlapped[[ncol + 4]], '_', 2)[, 1])  
  }
  dt <- merge(feature.dt, overlapped[, c('V1', 'V2', 'V3', 'V4', 'overlap', 'ph.id'), with = F], 
              by.x = c('chr', 'start', 'end', 'strand'), 
              by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  dt <- dt[dt$significant & dt$overlap, ]
  dt
}

CloseToExon <- function(tissue, feature, prefix){
  feature.path <- paste0(prefix, tissue, '_', feature, 'median_filtered_population.tsv')
  feature.dt <- fread(feature.path)
  x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
  system(x)  
  overlapped <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
  ncol <- dim(feature.dt)[2]
  overlapped$overlap <- T
  if(!any(grepl('ph.id', colnames(overlapped)))){
    overlapped$ph.id <- as.numeric(str_split_fixed(overlapped[[ncol + 4]], '_', 2)[, 1])  
  }
  dt <- merge(feature.dt, overlapped[, c('V1', 'V2', 'V3', 'V4', 'overlap', 'ph.id'), with = F], 
              by.x = c('chr', 'start', 'end', 'strand'), 
              by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  
  dt <- dt[dt$significant & dt$overlap, ]
  dt
}

SplicingScore <- function(ph){
  # RNA bridges
  exons.bridges <- fread(exons.bridges.path) 
  ph$RNAbridge <- F
  ph[ph$id %in% exons.bridges$ph_id, ]$RNAbridge <- T
  # Loop out DE exons
  dt <- LoopOut('kidney', 'exon', prefix)
  ph$kidney_loop_out_DE_exon <- F
  ph[ph$id %in% dt$ph.id, ]$kidney_loop_out_DE_exon <- T
  dt <- LoopOut('liver', 'exon', prefix)
  ph$liver_loop_out_DE_exon <- F
  ph[ph$id %in% dt$ph.id, ]$liver_loop_out_DE_exon <- T
  # Close to DE exon
  dt <- CloseToExon('kidney', 'exon', prefix)
  ph$kidney_close_to_DE_exon <- F
  ph[ph$id %in% dt$ph.id, ]$kidney_close_to_DE_exon <- T
  dt <- CloseToExon('liver', 'exon', prefix)
  ph$liver_close_to_DE_exon <- F
  ph[ph$id %in% dt$ph.id, ]$liver_close_to_DE_exon <- T
  ph$SplicingScore <- rowSums(ph[, c("RNAbridge", "kidney_loop_out_DE_exon", "liver_loop_out_DE_exon",  "kidney_close_to_DE_exon", "liver_close_to_DE_exon"), with = F])
  ph
}

ph <- SplicingScore(ph)
write.table(ph, path.to.ph, sep = '\t', row.names = F, quote = F)


