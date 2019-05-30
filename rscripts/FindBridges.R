library(data.table)
library(stringr)



MakeKDBed <- function(path.to.KD){
  KD.files <- list.files(path.to.KD, pattern = '*pooled', full.names = T)
  KD.list <- list()
  for (KD.file in KD.files){
    print(KD.file)
    KD.dt <- fread(KD.file)
    KD.dt <- KD.dt[abs(KD.dt$V3) > 0.1, ]
    KD.dt$chr <- str_split_fixed(KD.dt$V1, '_', 4)[, 1]
    KD.dt$start <- as.numeric(str_split_fixed(KD.dt$V1, '_', 4)[, 2])
    KD.dt$end <- as.numeric(str_split_fixed(KD.dt$V1, '_', 4)[, 3])
    KD.dt$strand <- str_split_fixed(KD.dt$V1, '_', 4)[, 4] 
    KD.dt$name <- paste(KD.dt$V2, '_dPSI=', KD.dt$V3, sep = '')
    KD.dt$score <- 1
    RBP <- unique(KD.dt$V2)
    KD.dt$thickStart <- KD.dt$start                                                                                                                                                                   
    KD.dt$thickSend <- KD.dt$end                                                                                                                                                                      
    KD.dt$itemRgb <- '255,0,0'
    KD.dt[KD.dt$V3 < 0, ]$itemRgb <- '0,0,255'
    KD.dt <- KD.dt[, c('chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickSend', 'itemRgb'), with = F]
    KD.list[[RBP]] <- KD.dt
  }
  KD.dt.all <- rbindlist(KD.list)
  write.table(KD.dt.all, '../RBPs/KD_all_RBP.bed', col.names = F, row.names = F, sep = '\t', quote = F)  
}

BedToBigBed <- function(path.to.bed, n_header_lines = 0){
  x <- paste0('tail -n +', n_header_lines + 1, ' ', path.to.bed, 
              ' | bedtools sort -i stdin > sorted.bed')
  system(x)
  x <- paste0('../tools/bedToBigBed sorted.bed ', '../tools/hg19.chrom.sizes ', gsub('bed', 'bb', path.to.bed))
  system(x)
  system('rm sorted.bed')
}

MakeFlankedHandles <- function(ph.handle.path.bed, flank.length = 500){
  handles <- fread(ph.handle.path.bed)
  handles$V2 <- handles$V2 - flank.length
  handles$V3 <- handles$V3 + flank.length
  write.table(handles, gsub('.bed', '_with_flanks.bed', ph.handle.path.bed), row.names = F, col.names = F, quote = F, sep = '\t')
}

IntersectKD <- function(path.to.handles, path.to.intervals, path.to.genes){
  x <- paste0('bedtools sort -i ', path.to.handles, 
              '| bedtools intersect -a stdin -b ', path.to.intervals, ' -wa -wb') 
  dt <- read.delim(pipe(x), header = F)
  dt$RBP <- str_split_fixed(dt$V8, '_', 2)[, 1]
  dt$ph_id <- as.numeric(str_split_fixed(dt$V4, '_', 2)[, 1])
  dt$RBP_handle <- str_split_fixed(dt$V4, '_', 2)[, 2]
  dt$RBP_dPSI <- as.numeric(str_split_fixed(dt$V8, '=', 2)[, 2])
  dt <- dt[, c('V1', 'V6', 'V7', 'ph_id', 'RBP_handle',  'RBP', 'RBP_dPSI')]
  write.table(dt, '../RBPs/ph_and_RBPs.tsv', sep = '\t', quote = F, row.names = F, col.names = F)
  x <- paste0('bedtools intersect -a ../RBPs/ph_and_RBPs.tsv -b ', path.to.genes, ' -wa -wb')
  dt <- read.delim(pipe(x), header = F)
  dt$gene_name <- gsub('gene_name ', '', dt$V12)
  dt <- dt[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'gene_name')]
  colnames(dt) <- c('chr', 'exon_start', 'exon_end', 'ph_id', 'RBP_handle', 'RBP', 'RBP_dPSI', 'gene_name')
  write.table(dt, '../RBPs/ph_and_RBPs.tsv', sep = '\t', quote = F, row.names = F)
}

MakeeClipBed <- function(path.to.eClip){
  dt <- fread(path.to.eClip)
  dt <- dt[dt$V7 > 3 & dt$V8 > 3,  ]
  dt$score <- 1
  dt <- dt[, c('V1', 'V2', 'V3', 'V4', 'score', 'V6'), with = F]
  write.table(dt, '../RBPs/eCLIP.bed', col.names = F, row.names = F, sep = '\t', quote = F)
}

IntersecteCLIP <- function(path.to.handles, path.to.intervals, path.to.genes){
  x <- paste0('bedtools sort -i ', path.to.handles, 
              '| bedtools intersect -a stdin -b ', path.to.intervals, ' -wa -wb') 
  dt <- read.delim(pipe(x), header = F)
  dt$RBP <- str_split_fixed(dt$V8, '_', 2)[, 1]
  dt$ph_id <- as.numeric(str_split_fixed(dt$V4, '_', 2)[, 1])
  dt$eCLIP_handle <- str_split_fixed(dt$V4, '_', 2)[, 2]
  dt <- dt[, c('V1', 'V6', 'V7', 'ph_id', 'eCLIP_handle', 'RBP')]
  write.table(dt, '../RBPs/ph_and_eCLIP.tsv', sep = '\t', quote = F, row.names = F, col.names = F)
  x <- paste0('bedtools intersect -a ../RBPs/ph_and_eCLIP.tsv -b ', path.to.genes, ' -wa -wb')
  dt <- read.delim(pipe(x), header = F)
  dt$gene_name <- gsub('gene_name ', '', dt$V11)
  dt <- dt[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'gene_name')]
  colnames(dt) <- c('chr', 'eCLIP_start', 'eCLIP_end', 'ph_id', 'eCLIP_handle', 'RBP', 'gene_name')
  write.table(dt, '../RBPs/ph_and_eCLIP.tsv', sep = '\t', quote = F, row.names = F)
}

IntersectKDeCLIP <- function(path.to.KD.ph, path.to.eClip.ph){
  KD <- fread(path.to.KD.ph)
  eCLIP <- fread(path.to.eClip.ph)
  merged <- merge(KD, eCLIP, by = c('chr', 'ph_id', 'RBP', 'gene_name'), all = F)
  merged <- merged[merged$RBP_handle != merged$eCLIP_handle, ]
  write.table(merged, '../RBPs/eCLIP_and_KD.tsv', sep = '\t', quote = F, row.names = F)
}

path.to.genes <- '../conservative_features/coding_genes.bed'
path.to.KD <- '../../dp/Projects/KD/shRNA/hg19/HepG2/'
MakeKDBed(path.to.KD)
path.to.KDbed <- '../RBPs/KD_all_RBP.bed' 
BedToBigBed(path.to.KDbed, n_header_lines = 0)
ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered_handles.bed'
MakeFlankedHandles(ph.handle.path.bed, flank.length = 500)
path.to.handles <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered_handles_with_flanks.bed'
IntersectKD(path.to.handles, path.to.KDbed, path.to.genes)
path.to.eClip <- '../../dp/ngs/encodedcc/hg19/eCLIP/joint_rep/HepG2.tsv'
MakeeClipBed(path.to.eClip)
path.to.eClip.bed <- '../RBPs/eCLIP.bed'
BedToBigBed(path.to.eClip.bed, n_header_lines = 0)
IntersecteCLIP(path.to.handles, path.to.eClip.bed, path.to.genes)
path.to.KD.ph <- "../RBPs/ph_and_RBPs.tsv"
path.to.eClip.ph <- "../RBPs/ph_and_eCLIP.tsv" 
IntersectKDeCLIP(path.to.KD.ph, path.to.eClip.ph)