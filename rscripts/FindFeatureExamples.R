library(data.table)



feature <- 'exon'
tissue <- 'kidney'


ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered_handles.bed'
ph.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered.bed'
ph.handle.flanked.path.bed <- gsub('.bed', '_with_flanks.bed', ph.handle.path.bed)
exons.bridges.path <- '../RBPs/eCLIP_and_KD.tsv'
prefix <- 'data/ipsa_output/ipsa_'
looped.out <- T

FindExamples <- function(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out, exons.bridges.path = ''){
  feature.path <- paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')
  print(feature.path)
  feature.dt <- fread(feature.path)
  if(feature == 'exon' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)  
  } else if(feature == 'exon' & looped.out == F){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)   
  } else if(feature == 'intron' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x) 
    dt <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
    ncol <- dim(feature.dt)[2]
    dt <- dt[dt$V2 <= dt[[ncol + 3]] & dt$V2 >= dt[[ncol + 2]] | dt$V3 <= dt[[ncol + 3]] & dt$V3 >= dt[[ncol + 2]], ]
    dt$ph.id <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 1]
    dt$handle <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 2]
    overlapped <- merge(dt, dt, by = c(paste("V", c(1:ncol), sep = ''), 'ph.id'))
    overlapped <- overlapped[overlapped$handle.y != overlapped$handle.x, ]
    overlapped$V21.x <- overlapped$V21.x + 1000
    overlapped$V21.y <- overlapped$V21.y + 1000
    overlapped$V22.x <- overlapped$V22.x - 1000
    overlapped$V22.y <- overlapped$V22.y - 1000
    overlapped <- overlapped[overlapped$V2 <= overlapped$V21.x & 
                               overlapped$V2 <= overlapped$V21.y & 
                               overlapped$V3 >= overlapped$V22.x & 
                               overlapped$V3 >= overlapped$V22.y,]
    to.drop <- grep('x|y', colnames(overlapped), value = T)
    overlapped <- overlapped[, ! names(overlapped) %in% c(to.drop), with = F]
    write.table(overlapped, paste0(prefix, tissue, '_', feature, '_overlapped.tsv'), row.names = F, col.names = T, quote = F, sep = '\t')
  }
  overlapped <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
  ncol <- dim(feature.dt)[2]
  overlapped$overlap <- T
  if(!any(grepl('ph.id', colnames(overlapped)))){
    overlapped$ph.id <- as.numeric(str_split_fixed(overlapped[[ncol + 4]], '_', 2)[, 1])  
  }
  dt <- merge(feature.dt, overlapped[, c('V1', 'V2', 'V3', 'V4', 'overlap', 'ph.id'), with = F], 
              by.x = c('chr', 'start', 'end', 'strand'), 
              by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  dt[is.na(dt$overlap), ]$overlap <- F
  
  if(grepl('site', feature)){
    dt <- dt[dt$end - dt$start <= 500, ]
  }
  dt$feature.id <- paste(dt$chr, dt$start, dt$end, sep = '_')
  
  if(feature == 'exon' & looped.out == F){
    exons.bridges <- fread(exons.bridges.path)
    exons.bridges$feature_id <- paste(exons.bridges$chr, exons.bridges$exon_start, exons.bridges$exon_end, sep = '_')
    dt$eCLIP <- F
    dt[dt$feature.id %in% exons.bridges$feature_id, ]$eCLIP <- T
  }
  
  dt <- dt[dt$significant & dt$overlap, ]
  write.table(dt, paste0(prefix, tissue, '_', feature, '_loopedOUT=', looped.out, '_examples.tsv'), sep = '\t', quote = F, row.names = F)
}

mtx <- data.frame(sign = c(length(unique(dt[dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[dt$significant & !dt$eCLIP, ]$feature.id))), 
                  n.sign = c(length(unique(dt[!dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[!dt$significant & !dt$eCLIP, ]$feature.id))))
print(mtx)
print(mtx$sign / (mtx$sign + mtx$n.sign) * 100)
fisher.test(mtx)$p.value

mtx <- data.frame(sign = c(length(unique(dt[dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[dt$significant & !dt$eCLIP & dt$overlap, ]$feature.id))), 
                  n.sign = c(length(unique(dt[!dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[!dt$significant & !dt$eCLIP &dt$overlap, ]$feature.id))))
print(mtx)
print(mtx$sign / (mtx$sign + mtx$n.sign) * 100)
fisher.test(mtx)$p.value




