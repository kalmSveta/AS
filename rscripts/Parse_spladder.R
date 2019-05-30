library(data.table)
library(stringr)
library(R.utils)

ParseSpladderOuput <- function(dt.path, anno.path, feature, tissue., out.path){
  dt <- fread(dt.path)
  anno <- fread(anno.path)
  
  if(feature == 'exon'){
    normal.columns <- grep(paste(anno[anno$tissue == tissue.,]$normal_file_id, collapse = '|'), colnames(dt), value = T)
    tumor.columns <- grep(paste(anno[anno$tissue == tissue.,]$tumor_file_id, collapse = '|'), colnames(dt), value = T)
    
    normal.inc.columns <- grep('intron_pre_conf|intron_aft_conf', normal.columns, value = T)
    tumor.inc.columns <- grep('intron_pre_conf|intron_aft_conf', tumor.columns, value = T)
    normal.exc.columns <- grep('intron_skip_conf', normal.columns, value = T)
    tumor.exc.columns <- grep('intron_skip_conf', tumor.columns, value = T)
    
    dt$inc_normal <- rowSums(dt[, normal.inc.columns, with = F])
    dt$inc_tumor <- rowSums(dt[, tumor.inc.columns, with = F])
    dt$exc_normal <- rowSums(dt[, normal.exc.columns, with = F])
    dt$exc_tumor <- rowSums(dt[, tumor.exc.columns, with = F])
    
    dt$chr <- paste('chr', dt$contig, sep = '')
    colnames(dt)[colnames(dt) == 'exon_start'] <- 'start'
    colnames(dt)[colnames(dt) == 'exon_end'] <- 'end'
    
    dt <- dt[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 'inc_normal'), with = F]
    write.table(dt, out.path, sep = '\t', row.names = F, quote = F)
  } else if(feature == 'intron'){
    normal.columns <- grep(paste(anno[anno$tissue == tissue.,]$normal_file_id, collapse = '|'), colnames(dt), value = T)
    tumor.columns <- grep(paste(anno[anno$tissue == tissue.,]$tumor_file_id, collapse = '|'), colnames(dt), value = T)
    
    normal.inc.columns <- grep('intron_cov', normal.columns, value = T)
    tumor.inc.columns <- grep('intron_cov', tumor.columns, value = T)
    normal.exc.columns <- grep('intron_conf', normal.columns, value = T)
    tumor.exc.columns <- grep('intron_conf', tumor.columns, value = T)
    
    dt$inc_normal <- rowSums(dt[, normal.inc.columns, with = F])
    dt$inc_tumor <- rowSums(dt[, tumor.inc.columns, with = F])
    dt$exc_normal <- rowSums(dt[, normal.exc.columns, with = F])
    dt$exc_tumor <- rowSums(dt[, tumor.exc.columns, with = F])
    
    dt$chr <- paste('chr', dt$contig, sep = '')
    colnames(dt)[colnames(dt) == 'intron_start'] <- 'start'
    colnames(dt)[colnames(dt) == 'intron_end'] <- 'end'
    
    dt <- dt[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 'inc_normal'), with = F]
    write.table(dt, out.path, sep = '\t', row.names = F, quote = F)
  } else if(feature == '3prime_site'){
    normal.columns <- grep(paste(anno[anno$tissue == tissue.,]$normal_file_id, collapse = '|'), colnames(dt), value = T)
    tumor.columns <- grep(paste(anno[anno$tissue == tissue.,]$tumor_file_id, collapse = '|'), colnames(dt), value = T)
    
    normal.inc.columns <- grep('intron1_conf', normal.columns, value = T)
    tumor.inc.columns <- grep('intron1_conf', tumor.columns, value = T)
    normal.exc.columns <- grep('intron2_conf', normal.columns, value = T)
    tumor.exc.columns <- grep('intron2_conf', tumor.columns, value = T)
    
    dt$inc_normal <- rowSums(dt[, normal.inc.columns, with = F])
    dt$inc_tumor <- rowSums(dt[, tumor.inc.columns, with = F])
    dt$exc_normal <- rowSums(dt[, normal.exc.columns, with = F])
    dt$exc_tumor <- rowSums(dt[, tumor.exc.columns, with = F])
    
    dt$chr <- paste('chr', dt$contig, sep = '')
    dt$start <- 0
    dt$end <- 0
    dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$start <- dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$exon_alt2_start
    dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$end <- dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$exon_alt1_start
    
    dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$start <- dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$exon_alt1_start
    dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$end <- dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$exon_alt2_start
    
    dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$start <- dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$exon_alt2_end
    dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$end <- dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$exon_alt1_end
    
    dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$start <- dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$exon_alt1_end
    dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$end <- dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$exon_alt2_end
    
    dt <- dt[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 'inc_normal'), with = F]
    write.table(dt, out.path, sep = '\t', row.names = F, quote = F)
  } else if(feature == '5prime_site'){
    normal.columns <- grep(paste(anno[anno$tissue == tissue.,]$normal_file_id, collapse = '|'), colnames(dt), value = T)
    tumor.columns <- grep(paste(anno[anno$tissue == tissue.,]$tumor_file_id, collapse = '|'), colnames(dt), value = T)
    
    normal.inc.columns <- grep('intron1_conf', normal.columns, value = T)
    tumor.inc.columns <- grep('intron1_conf', tumor.columns, value = T)
    normal.exc.columns <- grep('intron2_conf', normal.columns, value = T)
    tumor.exc.columns <- grep('intron2_conf', tumor.columns, value = T)
    
    dt$inc_normal <- rowSums(dt[, normal.inc.columns, with = F])
    dt$inc_tumor <- rowSums(dt[, tumor.inc.columns, with = F])
    dt$exc_normal <- rowSums(dt[, normal.exc.columns, with = F])
    dt$exc_tumor <- rowSums(dt[, tumor.exc.columns, with = F])
    
    dt$chr <- paste('chr', dt$contig, sep = '')
    dt$start <- 0
    dt$end <- 0
    dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$start <- dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$exon_alt2_start
    dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$end <- dt[dt$exon_alt1_start > dt$exon_alt2_start, ]$exon_alt1_start
    
    dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$start <- dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$exon_alt1_start
    dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$end <- dt[dt$exon_alt1_start < dt$exon_alt2_start, ]$exon_alt2_start
    
    dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$start <- dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$exon_alt2_end
    dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$end <- dt[dt$exon_alt1_end > dt$exon_alt2_end, ]$exon_alt1_end
    
    dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$start <- dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$exon_alt1_end
    dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$end <- dt[dt$exon_alt1_end < dt$exon_alt2_end, ]$exon_alt2_end
    
    dt <- dt[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 'inc_normal'), with = F]
    write.table(dt, out.path, sep = '\t', row.names = F, quote = F)
  }
}


features <- c('exon', 'intron', '5prime_site', '3prime_site')
features.names <- c('exon_skip', 'intron_retention', 'alt_5prime', 'alt_3prime')
tissues <- c('liver', 'kidney')
path <- 'data/spladder_output/'
anno.path <- '../delta_psi_icgc/metadata_normal-tumor.tsv'

for(tissue in tissues){
  for(i in c(1:length(features))){
    features.name <- features.names[i]
    feature <- features[i]
    dt.path <- paste0(path, '/merge_graphs_', tissue, '_', features.name, '_C3.txt.gz')   
    print(dt.path)
    out.path <- paste0(path, 'spladder_', tissue, '_', feature, '.tsv')
    ParseSpladderOuput(dt.path, anno.path, feature, tissue, out.path)
  }
}





