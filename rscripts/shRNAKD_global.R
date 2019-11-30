.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(reshape2)
library(stringr)
options(scipen = 999)

path.to.folder.anno <- '../../dp/ngs/encodedcc/'
anno <- paste0(path.to.folder.anno, 'shRNA_hg19.txt')
pairs <- paste0(path.to.folder.anno, 'shRNA_hg19_control_pairs.txt')
gene <- 'STAU1'
path.to.folder <- '/uge_mnt/home/dp/ipsa/hg19/shRNA/A07/'

MakeAnno <- function(anno, pairs, gene){
  dt <- fread(anno)
  KD <- subset(dt, grepl(gene, V13))
  KD <- KD[, c('V1', 'V4', 'V7', 'V13'), with = F]
  colnames(KD) <- c('KD.file', 'KD.patient', 'tissue', 'gene')
  pairs.dt <- fread(pairs, header = F)
  pairs.dt <- subset(pairs.dt, V1 %in% unique(KD$KD.patient))
  colnames(pairs.dt) <- c('KD.patient', 'CTR.patient')
  all <- merge(KD, pairs.dt, by = 'KD.patient')
  all <- merge(all, dt[, c('V1', 'V4'), with = F], by.x = 'CTR.patient', by.y = 'V4', all.x = T)
  colnames(all)[ncol(all)] <- 'CTR.file'
  all <- all[, !names(all) %in% c('CTR.patient', 'KD.patient'), with = F]
  all <- melt(all, id.vars = c('tissue', 'gene')) 
  all <- all[!duplicated(all),]
  all
}

ExtractIncExc <- function(path.to.folder, file){
  dt <- fread(paste0(path.to.folder, file, '.A07.gff'))
  dt <- subset(dt, V3 == 'exon')
  dt$inc <- as.numeric(gsub('inc', '', gsub('\\"', '', str_split_fixed(dt$V9, '; ', 5)[, 3])))
  dt$exc <- as.numeric(gsub('exc', '', gsub('\\"', '', str_split_fixed(dt$V9, '; ', 5)[, 2])))
  dt <- dt[, c('V1', 'V4', 'V5', 'V7', 'inc', 'exc'), with = F]
  colnames(dt) <- c('chr', 'start', 'stop', 'strand', paste0('inc_', file), paste0('exc_', file))
  dt
}


anno <- MakeAnno(anno, pairs, gene)

for(tissue. in unique(anno$tissue)){
  print(tissue.)
  anno.tissue <- subset(anno, tissue == tissue.)
  dt.list <- lapply(anno.tissue$value, function(file) ExtractIncExc(path.to.folder, file))
  dt <- Reduce(function(...) merge(..., by = c('chr', 'start', 'stop', 'strand')), dt.list)
  anno.tissue$tissue_var <- paste(anno.tissue$tissue, anno.tissue$var, sep = '_')
  # mean exc and inc, psi
  for(tissue.var in unique(anno.tissue$tissue_var)){
    files <- subset(anno.tissue, tissue_var == tissue.var)$value
    set(dt, j = paste0('inc_', tissue.var), 
        value = dt[, paste0('inc_', files[1]), with = F] +  dt[, paste0('inc_', files[2]), with = F])
    set(dt, j = paste0('exc_', tissue.var), 
        value = dt[, paste0('exc_', files[1]), with = F] +  dt[, paste0('exc_', files[2]), with = F])
    set(dt, j = paste0('psi_', tissue.var), 
        value = dt[, paste0('inc_', tissue.var), with = F] / 
          (dt[, paste0('inc_', tissue.var), with = F] + 2 * dt[, paste0('exc_', tissue.var), with = F]))
    dt <- subset(dt, get(paste0('inc_', tissue.var)) + 2 * get(paste0('exc_', tissue.var)) >= 20)
  }
  # delta psi
  set(dt, j = 'delta.psi', 
      value = dt[, paste0('psi_', tissue., '_KD.file'), with = F] - dt[, paste0('psi_', tissue., '_CTR.file'), with = F])
  # remove NA
  dt <- subset(dt, !is.na(delta.psi))
  # save
  dt$name <- c(1:nrow(dt))
  dt <- dt[, c('chr', 'start', 'stop', 'name', 'delta.psi', 'strand'), with = F]
  write.table(dt, paste0('../', gene, '/', 'delta_psi_', gene, '_', tissue., '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
}




