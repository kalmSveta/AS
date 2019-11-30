library(data.table)
library(stringr)
options(scipen = 999)

handles <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered_handles.bed'
MXE <- '../conservative_features/mut_ex_introns_human/introns_mut_ex_with_clusters.tsv'
ph <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered.tsv'

# if handle lies inside intron of cluster 
MarkMXE <- function(handles, MXE, ph){
  x <- paste0('bedtools sort -i ', handles, '| bedtools intersect -a stdin -b ', MXE, ' -u')
  overlap <- read.delim(pipe(x), header = F)
  overlap.ids <- str_split_fixed(overlap$V4, '_', 2)[, 1]
  ph <- fread(ph)
  ph$MXE <- FALSE
  ph[ph$id %in% overlap.ids, ]$MXE <- T
  ph
}

MXE <- '../conservative_features/mut_ex_introns_human/MXE.sg'
ph <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered.tsv'
ph.bed <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered.bed'

# If ph lies inside cluster
MarkMXE2 <- function(MXE, ph, ph.bed){
  # MXE <- fread(MXE)
  MXE$suppl <- 0
  MXE[MXE$V4 == '-', ]$suppl <- MXE[MXE$V4 == '-', ]$V2
  MXE[MXE$V4 == '-', ]$V2 <- MXE[MXE$V4 == '-', ]$V3
  MXE[MXE$V4 == '-', ]$V3 <- MXE[MXE$V4 == '-', ]$suppl
  MXE <- MXE[, ! names(MXE) %in% c('suppl'), with = F]
  
  write.table(MXE, '../conservative_features/mut_ex_introns_human/MXE_normal_order.tsv', sep = '\t', row.names = F, quote = F, col.names = F)
  MXE <- '../conservative_features/mut_ex_introns_human/MXE_normal_order.tsv'
  
  x <- paste0('bedtools sort -i ', ph.bed, '| bedtools intersect -a stdin -b ', MXE, ' -wa -wb')
  overlap <- read.delim(pipe(x), header = F)
  overlap <- overlap[overlap$V14 <= overlap$V2 & overlap$V3 <= overlap$V15, ]
  overlap.ids <- str_split_fixed(overlap$V4, '_', 2)[, 1]
  ph <- fread(ph)
  ph$MXE <- FALSE
  ph[ph$id %in% overlap.ids, ]$MXE <- T
  ph
}

ph.marked.MXE <- MarkMXE2(MXE, ph, ph.bed)
write.table(ph.marked.MXE, '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered.tsv', sep = '\t', row.names = F, quote = F)


dt <- read.table('../union_MXE_hg19.txt', sep = '\t', header = F, col.names = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6'), fill = T)