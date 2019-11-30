library(data.table)
library(stringr)
path <- '/uge_mnt/home/dp/ipsa/hg19/shRNA/A07/'
anno <- data.table(case = c(rep('KD', 4), rep('CTR', 4)), tissue = rep(c(rep('HepG2', 2), rep('K562', 2)), 2))
anno$replicate.id <- c('ENCFF423PGY', 'ENCFF666GBG', 'ENCFF396IPD', 'ENCFF340UTU', 'ENCFF236XCF', 'ENCFF901WHO', 'ENCFF109UAW', 'ENCFF827AVO')

dt.list <- list()
for(id in anno$replicate.id){
  dt <- fread(paste0(path, id, '.A07.gff'))
  dt <- dt[dt$V1=='chr6' & dt$V4>=36567344 & dt$V5 <= 36568807, ]
  dt$inc <- as.numeric(gsub('inc', '', gsub('\\"', '', str_split_fixed(dt$V9, '; ', 5)[, 3])))
  dt$exc <- as.numeric(gsub('exc', '', gsub('\\"', '', str_split_fixed(dt$V9, '; ', 5)[, 2])))
  dt$psi <- as.numeric(gsub('psi', '', gsub('\\"', '', str_split_fixed(dt$V9, '; ', 5)[, 4])))
  dt$denom <- dt$inc + 2 * dt$exc
  dt <- dt[, c('V1', 'V4', 'V5', 'inc', 'exc', 'psi', 'denom'), with = F]
  colnames(dt) <- c('chr', 'start', 'stop', paste0('inc_', id), paste0('exc_', id), paste0('psi_', id), paste0('denom_', id))
  dt.list[[id]] <- dt
}
dt <- Reduce(function(...) merge(..., by = c('chr', 'start', 'stop')), dt.list)


Stats <- function(tissue.){
  # KDs PSI
  ids <- anno[anno$tissue == tissue. & anno$case == 'KD', ]$replicate.id
  KDs.psi <- rowSums(dt[, paste0('psi_', ids), with = F]) / length(ids) # mean of two replicates
  KDs.denom <- rowSums(dt[, paste0('denom_', ids), with = F]) / length(ids) # mean of two replicates
  
  # CTRs PSI
  ids <- anno[anno$tissue == tissue. & anno$case == 'CTR', ]$replicate.id
  CTRs.psi <- rowSums(dt[, paste0('psi_', ids), with = F]) / length(ids) # mean of two replicates
  CTRs.denom <- rowSums(dt[, paste0('denom_', ids), with = F]) / length(ids) # mean of two replicates

  # delta PSI
  res <- dt[, c('chr', 'start', 'stop'), with = F]
  res$tissue <- tissue.
  res$dPSI <- KDs.psi - CTRs.psi
  res$ddenom <- KDs.denom - CTRs.denom
  res
}

res.K562 <- Stats('K562')
res.HepG2 <- Stats('HepG2')