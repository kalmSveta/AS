#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(stringr)
library(data.table)
options(scipen = 999)

# Usage: filter_ph.R path.to.ph path.to.ph.bed path.to.intervals n.ph.in.cluster. path.to.anno path.to.intervals path.to.tRNA path.to.snoRNA path.to.TFBS path.to.RM
args <- commandArgs(trailingOnly = TRUE)

# path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed.tsv'
# path.to.ph.bed <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed.bed'
# path.to.intervals <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/intervals_with_seqs.tsv'
n.ph.in.cluster. <- 15
path.to.anno <- '../conservative_features/gencode.v19.annotation.gtf'
path.to.tRNA <- '../conservative_features/tRNAs_hg19.txt'
path.to.snoRNA <- '../conservative_features/sno_miRNA_hg19.txt'
path.to.TFBS <- '../conservative_features/TFBS_hg19.txt'
path.to.RM <- '../conservative_features/RepeatMasker_hg19.txt'
path.to.CCDS <- '../conservative_features/ccds.txt'


path.to.ph <- args[1]
path.to.ph.bed <- args[2]
path.to.intervals <- args[3]
# n.ph.in.cluster. < args[4]
# path.to.anno <- args[5]
# path.to.tRNA <- args[6]
# path.to.snoRNA <- args[7]
# path.to.TFBS <- args[8]
# path.to.RM <- args[9]


SwapHandles <- function(path.to.ph.bed12){
  x <- paste('awk -F"\\t" \'{print $1,$11,$12,$4,$5,$6,$7,$8,$9,$10,$2,$3}\' OFS="\\t" ', 
             path.to.ph.bed12, 
             ' > ./swapt_handles.bed', sep = '')
  system(x)
  return(0)
}

IntersectTrack <- function(path.to.ph.bed12, path.to.track, track.name){
  if (track.name == 'TFBS'){
    min.overlap <- 0.5
  } else {
    min.overlap <- 1E-9
  }
  if (track.name == 'TFBS'){
    add <- paste('awk -F"\\t" \'$8 >= 2.33\' ', path.to.track, ' | ', sep = '')
    path.to.track <- ''
  } else{
    add <- ''
  }
  x <- paste(add, 'awk \'{printf ("%s\\t%s\\t%s\\n", $2, $3, $4)}\' ', 
             path.to.track, 
             '| tail -n +2| bedtools sort -i stdin | bedtools intersect -a ', 
             path.to.ph.bed12, 
             ' -b stdin -wa -wb -F ', min.overlap, sep = '')  
  dt <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  if(!is.null(dt)){
    return(dt$V4)  
  } else {
    return(NA)
  }
}

FilterTrack <- function(ph, path.to.ph.bed12, path.to.track, track.name){
  ids.left.handles.track <- unique(IntersectTrack(path.to.ph.bed12, path.to.track, track.name))
  SwapHandles(path.to.ph.bed12)
  ids.right.handles.track <- unique(IntersectTrack('swapt_handles.bed', path.to.track, track.name))
  ids.track <- intersect(ids.left.handles.track, ids.right.handles.track)
  set(ph, j = track.name, value = F)
  set(ph, i = which(ph$id %in% ids.track), j = track.name, value = T)
  print(paste('I have marked ph with ', track.name, '. Have killed ph:', sep = ''))
  print(dim(ph[which(ph[[track.name]] == T), ])[1])
  system('rm swapt_handles.bed')
  return(ph)
}

FilterCCDS <- function(ph, path.to.ph.bed12, path.to.track){
  x <- paste0('bedtools sort -i ', path.to.track, ' | bedtools intersect -a ', 
              path.to.ph.bed12, 
              ' -b stdin -wa -wb -f 1')   
  ids.left.handles.track <- unique(read.delim(pipe(x), header = F)$V4)
  SwapHandles(path.to.ph.bed12)
  x <- paste0('bedtools sort -i ', path.to.track, ' | bedtools intersect -a ', 
              'swapt_handles.bed', 
              ' -b stdin -wa -wb -f 1')  
  ids.right.handles.track <- unique(read.delim(pipe(x), header = F)$V4)
  ids.track <- unique(c(ids.left.handles.track, ids.right.handles.track))
  set(ph, j = 'CCDS', value = F)
  set(ph, i = which(ph$id %in% ids.track), j = 'CCDS', value = T)
  print(paste('I have marked ph with ', 'CCDS', '. Have killed ph:', sep = ''))
  print(dim(ph[which(ph[['CCDS']] == T), ])[1])
  system('rm swapt_handles.bed')
  return(ph)
}

IntersectRepeats <- function(path.to.ph.bed12, path.to.track){
  #add <- paste('awk -F"\\t" \'$12 == "Low_complexity" || $12 == "Simple_repeat"\' ', path.to.track, ' | ', sep = '')
  x <- paste('awk \'{printf ("%s\\t%s\\t%s\\t%s\\n", $6, $7, $8, $12)}\' ', 
             path.to.track,
             '| tail -n +2| bedtools sort -i stdin | bedtools intersect -a ', 
             path.to.ph.bed12, 
             ' -b stdin -wa -wb', sep = '') 
  print(cat(x))
  return(read.delim(pipe(x), header = F))
}

FilterRepeats <- function(ph, path.to.ph.bed12, path.to.track, track.name){
  df.left.handles.track <- IntersectRepeats(path.to.ph.bed12, path.to.track)
  SwapHandles(path.to.ph.bed12)
  df.right.handles.track <- IntersectRepeats('swapt_handles.bed', path.to.track)
  df.track <- rbind(df.left.handles.track[, c('V4', 'V16')], df.right.handles.track[, c('V4', 'V16')])
  df.track <- df.track[!duplicated(df.track$V4), ]
  ph <- merge(ph, df.track, by.x = 'id', by.y = 'V4', all.x = T)
  colnames(ph)[colnames(ph) == 'V16'] <- 'Repeat.type'
  ph[, (track.name) := F, ]
  set(ph, i = which(ph$Repeat.type %in% c('tRNA', 'Simple_repeat', 'rRNA', 'Low_complexity', 'snRNA','srpRNA','DNA')), j = track.name, value = T)
  print(paste('I have marked ph with ', track.name, '. Have killed ph:', sep = ''))
  x <- ph[[track.name]]
  print(dim(ph[x, ])[1])
  system('rm swapt_handles.bed')
  return(ph)    
}

IntersectCDS <- function(path.to.ph.bed12){
  x <- paste('bedtools sort -i ', 'anno_CDS_transcripts.bed', ' | bedtools intersect -a ', 
             path.to.ph.bed12, 
             ' -b stdin -wa -wb', sep = '')  
  return(read.delim(pipe(x), header = F)$V4)
}

FilterUntranslTranscripts <- function(ph, path.to.anno){
  anno <- fread(path.to.anno, header = F)
  anno.CDS <- anno[anno$V3 == 'CDS', ]
  CDS.transcript.id <- str_split_fixed(anno.CDS$V9, pattern = ';', Inf)[, 2]
  anno$transcript.id <- str_split_fixed(anno$V9, pattern = ';', Inf)[, 2] 
  anno <- anno[anno$transcript.id %in% CDS.transcript.id, ]
  anno <- anno[anno$V3 == 'transcript', ]
  write.table(anno[, c('V1', 'V4', 'V5', 'transcript.id')], 'anno_CDS_transcripts.bed', col.names = F, row.names = F, sep = '\t', quote = F)
  ids.CDS.left <- unique(IntersectCDS(path.to.ph.bed12))
  SwapHandles(path.to.ph.bed12)
  ids.CDS.right <- unique(IntersectCDS('swapt_handles.bed'))
  ph$CDS.transcript <- F
  ph[ph$id %in% ids.CDS.left & ph$id %in% ids.CDS.right, ]$CDS.transcript <- T
  print('I have marked ph with coding transcripts. Have killed ph:')
  print(dim(ph[ph$CDS.transcript == F, ])[1])
  system('rm swapt_handles.bed')
  system('rm anno_CDS_transcripts.bed')
  return(ph)
}



IntersectExons <- function(path.to.ph.bed12, path.to.anno, remove = 'exon', common.nts = 10){
  x <- paste0('tail -n +6 ', path.to.anno, 
             ' | awk -F"\\t" \'$3 == \"', remove, 
             '\"\'| awk \'{printf ("%s\\t%s\\t%s\\n", $1, $4, $5)}\' | bedtools sort -i stdin | bedtools intersect -a ',
             path.to.ph.bed12, ' -b stdin -wo ')
  dt <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  if(!is.null(dt)){
    dt <- dt[dt$V16 > common.nts, ]
    return(dt$V4)  
  } else {
    return(NA)
  }
}

FilterAlwaysIntronic <- function(ph, path.to.anno, remove = 'exon', common.nts = 10){
  ids.exons.left <- unique(IntersectExons(path.to.ph.bed12, path.to.anno, remove = remove, common.nts = common.nts))
  SwapHandles(path.to.ph.bed12)
  ids.exons.right <- unique(IntersectExons('swapt_handles.bed', path.to.anno, remove = remove, common.nts = common.nts))
  ph$always.intronic <- T
  ph[ph$id %in% ids.exons.left | ph$id %in% ids.exons.right, ]$always.intronic <- F
  print('I have marked ph with always intronic. Have killed ph:')
  print(dim(ph[ph$always.intronic == F, ])[1])
  system('rm swapt_handles.bed')
  return(ph)
}

FilterSubopt <- function(ph, path.to.intervals){
  x <- paste0('tail -n +2 ', path.to.intervals, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', path.to.ph.bed12, ' -wa -wb')
  dt <- read.delim(pipe(x), header = F)
  left <- dt[, c('V8', 'V12')]
  colnames(left) <- c('left.interval.id', 'ph.id')
  SwapHandles(path.to.ph.bed12)
  x <- paste0('tail -n +2 ', path.to.intervals, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', 'swapt_handles.bed', ' -wa -wb')
  dt <- read.delim(pipe(x), header = F)
  
  right <- dt[, c('V8', 'V12')]
  colnames(right) <- c('right.interval.id', 'ph.id')
  
  dt <- merge(ph, left, by.x = 'id', by.y = 'ph.id')
  dt <- merge(dt, right, by.x = 'id', by.y = 'ph.id')
  
  dt <- dt[order(dt$energy), ]
  dt <- dt[! duplicated(dt[, c('left.interval.id', 'right.interval.id'), with = F]), ]
  ph$suboptimal <- T
  ph[ph$id %in% dt$id, ]$suboptimal <- F
  system('rm swapt_handles.bed')
  return(ph)
}

Make_bed12 <- function(path.to.ph){
  ph <- fread(path.to.ph)
  ph_bed = ph[, c("chr", "panhandle_start", "panhandle_left_hand", "id", "energy", "strand", "alignment1",
               "alignment2", "al1_length", "al2_length", "panhandle_right_hand", "panhandle_end"), with = F]
  write.table(ph_bed, gsub('.tsv', '.bed12', path.to.ph), row.names = F, col.names = F, sep = "\t", quote = F)
}

Count_clusters <- function(ph, path.to.ph.bed){
  system(paste0('tail -n +2 ', path.to.ph.bed, '>./tmp'))
  ph.tmp <- read.delim(pipe('bedtools intersect -a ./tmp -b ./tmp -c'), header = F)
  system('rm ./tmp')
  ph.tmp$id <- as.numeric(str_split_fixed(ph.tmp$V4, '_', 2)[, 1])
  ph <- merge(ph, ph.tmp[, c('id', 'V13')], by = 'id')
  colnames(ph)[colnames(ph) == 'V13'] <- 'n.ph.in.cluster' 
  ph
}

MakePlusStrandSeq <- function(ph){
  ph$alignment1_C <- unlist(lapply(ph$alignment1, function(x) chartr("ATGC","TACG", x)))
  ph$alignment2_C <- unlist(lapply(ph$alignment2, function(x) chartr("ATGC","TACG", x)))
  ph[ph$strand == '-', ]$alignment1 <- ph[ph$strand == '-', ]$alignment1_C
  ph[ph$strand == '-', ]$alignment2 <- ph[ph$strand == '-', ]$alignment2_C
  ph <- ph[, !names(ph) %in% c('alignment1_C', 'alignment2_C'), with = F]
  ph
}

Make_bed12(path.to.ph)
print('Made bed12')
path.to.ph.bed12 <- gsub('.tsv', '.bed12', path.to.ph)
ph <- fread(path.to.ph)
print('read in file')
ph <- FilterRepeats(ph, path.to.ph.bed12, path.to.RM, 'Repeats')
ph <- FilterTrack(ph, path.to.ph.bed12, path.to.tRNA, 'tRNA')
ph <- FilterTrack(ph, path.to.ph.bed12, path.to.snoRNA, 'sno.miRNA')
ph <- FilterTrack(ph, path.to.ph.bed12, path.to.TFBS, 'TFBS')
ph <- FilterUntranslTranscripts(ph, path.to.anno)
#ph <- FilterAlwaysIntronic(ph, path.to.anno, remove = 'CDS', common.nts = 10)
#ph <- MakePlusStrandSeq(ph)
#ph <- FilterCCDS(ph, path.to.ph.bed12, path.to.CCDS)
#ph <- FilterSubopt(ph, path.to.intervals)
#ph <- Count_clusters(ph, path.to.ph.bed)

filtered <- ph[ # ph$suboptimal == F & 
                # ph$CCDS == T &
                # ph$always.intronic == T & 
                # ph$CDS.transcript == T &
                 ph$tRNA == F & 
                 ph$sno.miRNA == F & 
                 ph$TFBS == F & 
                 ph$Repeats == F, ]
print('I have ph after filtering:')
print(dim(filtered)[1])

write.table(ph, paste(path.to.ph, '2', sep = ''), sep = '\t', row.names = F, quote = F)
write.table(filtered, gsub('.tsv', '_filtered.tsv', path.to.ph), sep = '\t', row.names = F, quote = F)

# filtered2 <- filtered[filtered$n.ph.in.cluster < n.ph.in.cluster., ]
# write.table(filtered2, gsub('.tsv', '_filtered2.tsv', path.to.ph), sep = '\t', row.names = F, quote = F)
# print('I have ph in small clusters:')
# print(dim(filtered2)[1])

