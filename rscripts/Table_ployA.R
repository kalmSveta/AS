library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)
library(Cairo)
library(RColorBrewer)

SelectWorkingSample <- function(dt, genes){
  new.dt.name <- gsub('.bed', '_sample.bed', dt)
  system(paste0('bedtools intersect -a ', dt, ' -b ', genes, ' -f 0.9 -u -s > ', new.dt.name))
  return(new.dt.name)
}


genes <- '../conservative_features/not_intersected_coding_genes.bed'

#introns <- '../introns/introns_K562.bed'
#introns <- '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed'
introns <- '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_balanced_set_for_tables.bed'

introns <- SelectWorkingSample(introns, genes)

path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')


# ends
polyA <- '../polyA/GSE30198_human.pas.uniq_filtered.bed'
transcript.ends <- '../conservative_features/gencode.v19.annotation_transcripts_ends.bed'

# starts
CAGE <- '../CAGE/NHEK_polyA+_all.bed'
transcript.starts <- '../conservative_features/gencode.v19.annotation_transcripts_starts.bed'

energy.cutoffs <- c(-15, -20, -25, -30)
energy.colors <- brewer.pal(n = 12, "Paired")[c(4, 7, 8, 6)]
names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')


Count <- function(introns, intervals1.path, intervals2.path, intervals1.name, intervals2.name, energy.cutoff = F){
  dt.intervals1 <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', introns, ' -b ', intervals1.path, ' -s -wa -wb -F 1')), header = F))
  dt.intervals2 <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', introns, ' -b ', intervals2.path, ' -s -wa -wb -F 1')), header = F))
  dt <- fread(introns)
  
  if(energy.cutoff){
    dt.intervals1$energy <- as.numeric(str_split_fixed(dt.intervals1$V10, '_', 2)[, 2])
    dt.intervals1 <- subset(dt.intervals1, energy <= energy.cutoff)
  }
  dt$id <- paste(dt$V1, dt$V2, dt$V3, dt$V6, sep = '_')
  dt.intervals2$id <- paste(dt.intervals2$V1, dt.intervals2$V2, dt.intervals2$V3, dt.intervals2$V6, sep = '_')
  dt.intervals1$id <- paste(dt.intervals1$V1, dt.intervals1$V2, dt.intervals1$V3, dt.intervals1$V6, sep = '_')
  dt$interval1 <- F
  dt[dt$id %in% dt.intervals1$id, ]$interval1 <- T
  dt$interval2 <- F
  dt[dt$id %in% dt.intervals2$id, ]$interval2 <- T
  
  agr <- aggregate(V5 ~ interval1 + interval2, dt, FUN = sum)

  result <- dcast(agr, interval1 ~ interval2)
  result$interval1 <- gsub('TRUE', intervals1.name, result$interval1)
  result$interval1 <- gsub('FALSE', paste0('no ', intervals1.name), result$interval1)
  rownames(result) <- result$interval1
  result <- result[, c('FALSE', 'TRUE')]
  colnames(result) <- c(paste0('no ', intervals2.name), intervals2.name)
  
  result
}

Percentage <- function(dt){
  perc <- apply(dt, 1, function(row) 
    paste(round(as.numeric(row[2]) / (as.numeric(row[1]) + as.numeric(row[2]))*100, 2), '%', sep = ''))
  dt[, 2] <- paste(dt[, 2], '(', perc, ')', sep = '')
  dt
}

Cond_prob <- function(introns, phs, ends, starts, energy.cutoff = -15){
  dt.phs <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', introns, ' -b ', phs, ' -s -wa -wb -F 1')), header = F))
  dt.ends <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', introns, ' -b ', ends, ' -s -wa -wb -F 1')), header = F))
  dt.starts <- as.data.table(read.delim(pipe(paste0('bedtools intersect -a ', introns, ' -b ', starts, ' -s -wa -wb -F 1')), header = F))
  dt <- fread(introns)
  
  dt.phs$energy <- as.numeric(str_split_fixed(dt.phs$V10, '_', 2)[, 2])
  dt.phs <- subset(dt.phs, energy <= energy.cutoff) 
  
  dt$id <- paste(dt$V1, dt$V2, dt$V3, dt$V6, sep = '_')
  dt.phs$id <- paste(dt.phs$V1, dt.phs$V2, dt.phs$V3, dt.phs$V6, sep = '_')
  dt.starts$id <- paste(dt.starts$V1, dt.starts$V2, dt.starts$V3, dt.starts$V6, sep = '_')
  dt.ends$id <- paste(dt.ends$V1, dt.ends$V2, dt.ends$V3, dt.ends$V6, sep = '_')
  
  dt$ph <- F
  dt[dt$id %in% dt.phs$id, ]$ph <- T
  dt$start <- F
  dt[dt$id %in% dt.starts$id, ]$start <- T
  dt$end <- F
  dt[dt$id %in% dt.ends$id, ]$end <- T
  
  result <- lapply(c(FALSE, TRUE), function(ph.){
    tmp <- subset(dt, ph == ph.)
    agr <- aggregate(V5 ~ start + end, tmp, FUN = sum)
    c(paste(round(subset(agr, start == T & end == F)$V5 / sum(subset(agr, end == F)$V5) * 100, 2), '%', sep = ''),
    paste(round(subset(agr, start == T & end == T)$V5 / sum(subset(agr, end == T)$V5) * 100, 2), '%', sep = ''))
  })
  result <- data.frame(result)
  colnames(result) <- c('ph-', 'ph+')
  rownames(result) <- c('start | no end', 'start | end')
  result
}

cairo_pdf('phs_introns_transcripts_simple_introns_anno.pdf', onefile = T)
for(energy.cutoff. in energy.cutoffs){
  phs.vs.ends <- Count(introns, intervals1.path = phs, intervals2.path = transcript.ends, 
                       intervals1.name = 'ph', intervals2.name = 'end', energy.cutoff = energy.cutoff.)
  phs.vs.ends.pretty <- Percentage(phs.vs.ends)
  
  phs.vs.starts <- Count(introns, intervals1.path = phs, intervals2.path = transcript.starts, 
                         intervals1.name = 'ph', intervals2.name = 'start', energy.cutoff = energy.cutoff.)
  phs.vs.starts.pretty <- Percentage(phs.vs.starts)
  
  ends.vs.starts <- Count(introns, intervals1.path = transcript.ends, intervals2.path = transcript.starts, 
                          intervals1.name = 'end', intervals2.name = 'start')
  ends.vs.starts.pretty <- Percentage(ends.vs.starts)
  
  all <- Cond_prob(introns, phs, transcript.ends, transcript.starts, energy.cutoff = energy.cutoff.)  
  
  grid.arrange(
    tableGrob(phs.vs.ends.pretty),
    tableGrob(phs.vs.starts.pretty),
    tableGrob(ends.vs.starts.pretty),
    tableGrob(all),
    nrow=2, ncol = 2, top = paste0('\u0394G<=', energy.cutoff., 'kcal/mol'))
}
dev.off()

