library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(qvalue)
library(stringr)
library(data.table)
source('functions.R')
options(scipen = 999)


Preprocess <- function(dt){
  colnames(dt) <- gsub('sum_', '', colnames(dt))
  colnames(dt) <- gsub('median_', '', colnames(dt))
  colnames(dt) <- gsub('ret', 'exc', colnames(dt))
  colnames(dt) <- gsub('teta', 'psi', colnames(dt))
  colnames(dt) <- gsub('stop', 'end', colnames(dt))
  if (any(grepl('coord', colnames(dt)))){
    dt$start <- dt$coord
    dt$end <- dt$coord + 1
  } else if(!any(colnames(dt) == 'chr')){
    dt$chr <- str_split_fixed(dt$exon_chr_start_end_strand, ':', Inf)[, 1]
    dt$strand <- unlist(lapply(str_split_fixed(str_split_fixed(dt$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 2], function(x) substr(x, nchar(x), nchar(x))))
    dt$start <- as.numeric(str_split_fixed(str_split_fixed(dt$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 1])
    dt$end <- as.numeric(unlist(lapply(str_split_fixed(str_split_fixed(dt$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 2], function(x) substr(x, 1, nchar(x) - 1))))
  } 
  dt$exc_tumor <- as.numeric(dt$exc_tumor)
  dt$inc_tumor <- as.numeric(dt$inc_tumor)
  dt$exc_normal <- as.numeric(dt$exc_normal)
  dt$inc_normal <- as.numeric(dt$inc_normal)
  dt <- dt[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 
               'inc_normal'), with = F]
  dt
}

Aggregate <- function(dt){
  dt.new <- aggregate(dt[, c('exc_tumor', 'inc_tumor', 'exc_normal', 'inc_normal'), with = F], by = list(dt$chr, dt$start, dt$end, dt$strand), FUN=sum)
  colnames(dt.new) <- colnames(dt)
  dt.new
}

SelectCodingExons <- function(dt, path.to.CDS.exons){
  write.table(dt, 'tmp.bed', sep = '\t', col.names = F, row.names = F, quote = F)
  columns <- colnames(dt)
  x <- paste0('bedtools sort -i tmp.bed | bedtools intersect -f 1 -r -a stdin -b ', path.to.CDS.exons, ' -wa -wb')
  dt <- as.data.table(read.delim(pipe(x), header = F))
  dt <- dt[, 1:(ncol(dt) - 4)]
  colnames(dt) <- columns
  dt <- dt[! duplicated(dt), ]
  x <- 'rm tmp.bed'
  system(x)
  dt
}

CalculateSJcount <- function(dt, feature){ 
  if(feature == 'exon'){ # c(inc, exc)
    coef <- c(1, 2)
  } else if(grepl('site', feature)){ # c(inc, ret)
    coef <- c(1, 1)
  } else if(feature == 'intron'){ # c(inc, ret)
    coef <- c(2, 1)
  }
  dt$sjcount_normal <- coef[1] * dt$inc_normal + coef[2] * dt$exc_normal 
  dt$sjcount_tumor <- coef[1] * dt$inc_tumor + coef[2] * dt$exc_tumor
  dt
}

CalculatePsi <- function(dt){
  dt$psi_normal <- dt$inc_normal / dt$sjcount_normal
  dt$psi_tumor <- dt$inc_tumor / dt$sjcount_tumor
  dt$delta_psi <- dt$psi_tumor - dt$psi_normal
  dt
}
    
SelectDetectableFeatures <- function(dt, denom.cutoff){
  dt <- dt[dt$sjcount_normal >= denom.cutoff & dt$sjcount_tumor >= denom.cutoff, ]
  dt  
}

RegressOutExpression <- function(dt){
  dt$delta_psi = round(dt$delta_psi, digits = 2)
  dt$log10sjcount = round(with(dt, log10(sjcount_tumor + sjcount_normal)), digits = 4)
  dt$log10FC = round(with(dt, log10(sjcount_tumor / sjcount_normal)), digits = 4)
  model = lm(delta_psi ~ log10FC, dt)
  #print(summary(model))
  dt$delta_psi_c = round(model$residuals, digits = 2) + p.cutoff / 2
  dt
}

IntersectWithGenes <- function(dt, path.to.coding.genes){
  write.table(dt, 'tmp.bed', sep = '\t', col.names = F, row.names = F, quote = F)
  genes <- read.delim(pipe(paste("bedtools intersect -a tmp.bed -b ", path.to.coding.genes, " -wa -wb", sep = '')), header = F)
  genes <- genes[genes$V4 == genes$V22, ]
  genes <- genes[!duplicated(genes[, c("V1", 'V2', "V3"), ]), ]
  genes$V21 <- gsub('gene_name ', '', genes$V21)
  dt <- merge(dt, genes[, c('V1', 'V2', 'V3', 'V4', 'V21')], by.x = c('chr', 'start', 'end', 'strand'), by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  colnames(dt)[colnames(dt) == 'V21'] <- 'gene.name'
  x <- 'rm tmp.bed'
  system(x)
  dt <- dt[!is.na(dt$gene.name), ]
  dt
}

CalculateQval <- function(dt, q.cutoff, p.cutoff, do.plot = F){
  dt1 = subset(dt, delta_psi != 0)
  dt1$bin = cut(dt1$log10sjcount, breaks = seq(0, round(max(dt1$log10sjcount), digits = 0), 0.5))
  merge(dt1, ddply(dt1,.(bin), summarise, mean = mean(delta_psi_c),
                   sd = sd(delta_psi_c)), by = 'bin') -> dt2
  dt2$z = with(dt2, round((delta_psi_c - mean) / sd, digits = 2))
  dt2$p = pnorm(-abs(dt2$z))
  dt2$q = qvalue_truncp(dt2$p)$qvalues
  
  if(do.plot){
    dt3 = subset(dt2, q < q.cutoff) %>% group_by(gene.name) %>% slice(which.min(q))
    p = ggplot(dt3, aes(x = log10sjcount, y = delta_psi_c)) +
      geom_point(size = 0.5, alpha = 0.5, aes(color = (q < q.cutoff))) + 
      geom_abline(slope = 0,lty = "dashed") + 
      geom_label_repel(size = 3, data = dt3, aes(label = gene.name, color = (q < q.cutoff)), 
                       min.segment.length = 0, force = 10) + 
      xlab("log(SJ count)") +
      ylab(expression(paste(Delta, Psi,"'(Tumor-Normal)"))) + 
      theme_classic() + theme(legend.position = "none") +
      scale_color_brewer(palette = "Set2") + 
      ggtitle(title)
    
    png(gsub('.tsv', '_qvalue.png', input), width = 6, height = 7, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  dt2 <- as.data.table(dt2)
  dt.new <- merge(dt, dt2[, c("chr", "start", "end", "strand", "q"), with = F], all.x = T, by = c("chr", "start", "end", "strand"))
  dt.new$significant <- F
  dt.new[!is.na(dt.new$q) & dt.new$q < q.cutoff, ]$significant <- T
  dt.new
}

SaveTrack <- function(dt, tissue, feature, prefix, title = ''){
  to.save <- as.data.table(dt[dt$significant, ])
  to.save$index <- paste('deltaPSI=', to.save$delta_psi_c, '|inc_normal=', to.save$inc_normal, '|exc_normal=', to.save$exc_normal, 
                                                         '|inc_tumor=', to.save$inc_tumor, '|exc_tumor=', to.save$exc_tumor, sep = '')
  to.save$score <- 0
  to.save$color <- '0,0,255'
  to.save[to.save$delta_psi > 0, ]$color <- '255,0,0'
  to.save$thick_start <- to.save$start
  to.save$thick_end <- to.save$end
  to.save <- to.save[, c('chr', 'start', 'end', 'index', 'score', 'strand', 'thick_start', 'thick_end', 'color'), with = F]
  path.to.bed <- paste0(prefix, tissue, '_', feature, '.bed')
  write.table(to.save, path.to.bed, sep = '\t', row.names = F, col.names = F, quote = F)
  BedToBigBed(path.to.bed, n_header_lines = 0)
} 

main <- function(feature, tissue, path.to.coding.genes, path.to.CDS.exons, q.cutoff, denom.cutoff, title, prefix){
  input <- paste0(prefix, tissue, '_', feature, '.tsv')
  dt <- fread(input)
  dt <- Preprocess(dt)
  dt <- Aggregate(dt)
  if(feature == 'exon'){
    dt <- SelectCodingExons(dt, path.to.CDS.exons)
  }
  dt <- CalculateSJcount(dt, feature)
  dt <- CalculatePsi(dt)
  dt <- SelectDetectableFeatures(dt, denom.cutoff)
  dt <- RegressOutExpression(dt)
  dt <- IntersectWithGenes(dt, path.to.coding.genes)
  dt <- CalculateQval(dt, q.cutoff, p.cutoff, do.plot = F)
  write.table(dt, paste0(prefix, tissue, '_', feature, '_filtered_population.tsv'), sep = '\t', row.names = F, quote = F) 
  SaveTrack(dt, tissue, feature, prefix)
}  

prefix <- 'data/spladder_output/spladder_'
denom.cutoff <- 30
q.cutoff <- 0.05
p.cutoff <- 0.05
path.to.CDS.exons <- '../conservative_features/coding_exons.bed'
path.to.coding.genes <-  "../conservative_features/coding_genes.bed"
features <- c('5prime_site', '3prime_site', 'exon', 'intron')
#features <- c('site', 'exon', 'intron')

for (tissue in c('kidney', 'liver')){
  for (feature in features){
    title <- paste0(feature, '_', tissue)
    print(title)
    main(feature, tissue, path.to.coding.genes, path.to.CDS.exons, q.cutoff, denom.cutoff, title, prefix)
  }
}

prefix <- 'data/ipsa_output/ipsa_'
for (tissue in c('kidney', 'liver')){
  for (feature in c('site', 'exon', 'intron')){
    title <- paste0(feature, '_', tissue)
    print(title)
    dt <- fread(paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')) 
    SaveTrack(dt, tissue, feature, prefix)
  }
}
prefix <- 'data/spladder_output/spladder_'
for (tissue in c('kidney', 'liver')){
  for (feature in c('5prime_site', '3prime_site', 'exon', 'intron')){
    title <- paste0(feature, '_', tissue)
    print(title)
    dt <- fread(paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')) 
    SaveTrack(dt, tissue, feature, prefix)
  }
}


