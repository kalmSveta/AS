library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(qvalue)
library(stringr)
library(data.table)
source('functions.R')
options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

# Usage: Rscript qvalue.R input title q.cutoff denom.cutoff feature


input <- 'data/kidney_exon.tsv'
feature <- 'exon'
stat <- 'qvalue'
loop.out <- T
path.to.chrom.file <- '../tools/hg19.chrom.sizes'
path.to.ph.bed12 <- '../conservative_features/whole_human_genome/alignments_whole_human_genome_processed_bed12.bed'
path.to.ph.main <- '../conservative_features/whole_human_genome/alignments_whole_human_genome_processed.tsv'
path.to.mut <- '../mutations_icgc/SNVs/RECA_simple_somatic_mutation.open_only_SNV.bed'
path.to.CDS.exons <- '../conservative_features/coding.exons.bed'
q.cutoff <- 0.05
p.cutoff <- 0.05
delta.psi.cutoff <- 0.07
denom.cutoff <- 30
flank.length <- 100



Preprocess <- function(df){
  colnames(df) <- gsub('sum_', '', colnames(df))
  colnames(df) <- gsub('ret', 'exc', colnames(df))
  colnames(df) <- gsub('teta', 'psi', colnames(df))
  colnames(df) <- gsub('stop', 'end', colnames(df))
  if (feature == 'site'){
    df <- df[, c('chr', 'coord', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 
                 'inc_normal'), with = F]
  } else if(feature == 'exon'){
    df$chr <- str_split_fixed(df$exon_chr_start_end_strand, ':', Inf)[, 1]
    df$start <- str_split_fixed(str_split_fixed(df$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 1]
    df$end <- unlist(lapply(str_split_fixed(str_split_fixed(df$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 2], function(x) substr(x, 1, nchar(x) - 1)))
    df$strand <- unlist(lapply(str_split_fixed(str_split_fixed(df$exon_chr_start_end_strand, ':', Inf)[, 2], '_', Inf)[, 2], function(x) substr(x, nchar(x), nchar(x))))
    df$start <- as.numeric(df$start)
    df$end <- as.numeric(df$end)
    df <- df[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 
                 'inc_normal'), with = F]
  } else if(feature == 'intron'){
    df <- df[, c('chr', 'start', 'end', 'strand', 'exc_tumor', 'inc_tumor', 'exc_normal', 
                 'inc_normal'), with = F]
  } else{
    print('Do not know such feature')
  }
  df
}

CalculateSJcount <- function(df, feature){ 
  if(feature == 'exon'){ # c(inc, exc)
    coef <- c(1, 2)
  } else if(feature == 'site'){ # c(inc, ret)
    coef <- c(1, 1)
  } else if(feature == 'intron'){ # c(inc, ret)
    coef <- c(2, 1)
  }
  df$sjcount_normal <- coef[1] * df$inc_normal + coef[2] * df$exc_normal 
  df$sjcount_tumor <- coef[1] * df$inc_tumor + coef[2] * df$exc_tumor
  df
}

CalculatePsi <- function(df){
  df$psi_normal <- df$inc_normal / df$sjcount_normal
  df$psi_tumor <- df$inc_tumor / df$sjcount_tumor
  df$delta_psi <- df$psi_tumor - df$psi_normal
  df
}

IntersectWithGenes <- function(df, denom.cutoff){
  df.bed <- df
  df.bed$name <- '.'
  df.bed$score <- '.'
  if (feature == 'site'){
    df.bed <- df.bed[, c('chr', 'coord', 'coord', 'name', 'score', 'strand'), with = F]   
  } else {
    df.bed <- df.bed[, c('chr', 'start', 'end', 'name', 'score', 'strand'), with = F]   
  }
  write.table(df.bed, gsub('.tsv', paste('_filtered_denominator', denom.cutoff, '.bed', sep = ''), input), sep = '\t', row.names = F, col.names = F, quote = F)
  genes <- read.delim(pipe(paste("bedtools intersect -a ", 
                                 gsub('.tsv', paste('_filtered_denominator', 
                                                    denom.cutoff, '.bed', sep = ''), input), 
                                 " -b ../conservative_features/gencode.v19.annotation_protein_coding_genes_gene_names.bed -wa -wb", sep = '')), 
                      header = F)
  if(feature == 'site'){
    df <- merge(df, genes[, c('V1', 'V2', 'V6', 'V10')], by.x = c('chr', 'strand', 'coord'), by.y = c('V1', 'V6', 'V2'), all = F)
    df <- df[!duplicated(df[, c('chr', 'strand', 'coord'), with = F]), ]
  } else {
    df <- merge(df, genes[, c('V1', 'V2', 'V3', 'V6', 'V10')], by.x = c('chr', 'strand', 'start', 'end'), by.y = c('V1', 'V6', 'V2', 'V3'), all = F)
    df <- df[!duplicated(df[, c('chr', 'strand', 'start', 'end'), with = F]), ]
  }
  colnames(df)[colnames(df) == 'V10'] <- 'gene.name'
  df
}

RegressOutExpression <- function(df){
  df$delta_psi = round(df$delta_psi, digits = 2)
  df$log10sjcount = round(with(df, log10(sjcount_tumor + sjcount_normal)), digits = 4)
  df$log10FC = round(with(df, log10(sjcount_tumor / sjcount_normal)), digits = 4)
  model = lm(delta_psi ~ log10FC, df)
  #print(summary(model))
  
  df$delta_psi_c = round(model$residuals, digits = 2) + p.cutoff / 2
  df
}

CalculateQval <- function(df, q.cutoff, p.cutoff, denom.cutoff){
  df1 = subset(df, delta_psi != 0)
  df1$bin = cut(df1$log10sjcount, breaks = seq(0, round(max(df1$log10sjcount), digits = 0), 0.5))
  merge(df1, ddply(df1,.(bin), summarise, mean = mean(delta_psi_c),
                   sd = sd(delta_psi_c)), by = 'bin') -> df2
  
  df2$z = with(df2, round((delta_psi_c - mean) / sd, digits = 2))
  df2$p = pnorm(-abs(df2$z))
  df2$q = qvalue_truncp(df2$p)$qvalues
  df3 = subset(df2, q < q.cutoff) %>% group_by(gene.name) %>% slice(which.min(q))
  
  p = ggplot(df2, aes(x = log10sjcount, y = delta_psi_c)) +
    geom_point(size = 0.5, alpha = 0.5, aes(color = (q < q.cutoff))) + 
    geom_abline(slope = 0,lty = "dashed") + 
    geom_label_repel(size = 3, data = df3, aes(label = gene.name, color = (q < q.cutoff)), 
                     min.segment.length = 0, force = 10) + 
    xlab("log(SJ count)") +
    ylab(expression(paste(Delta, Psi,"'(Tumor-Normal)"))) + 
    theme_classic() + theme(legend.position = "none") +
    scale_color_brewer(palette = "Set2") + 
    ggtitle(title)
  
  png(gsub('.tsv', '_qvalue.png', input), width = 6, height = 7, units = 'in', res = 300)
  print(p)
  dev.off()
  
  write.table(df2, gsub('.tsv', paste('filtered_denominator', denom.cutoff, '_qval.tsv', sep = ''), input), sep = '\t', row.names = F, quote = F)
  df2
}

SwapHandles <- function(path.to.ph.bed12){
  x <- paste('awk -F"\\t" \'{print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}\' OFS="\\t" ', 
             path.to.ph.bed12, 
             ' > ./swapt_handles.bed', sep = '')
  system(x)
  return(0)
}

IntersectHandle <- function(path.to.ph.bed12, path.to.chrom.file, path.to.feature.stat, feature, flank.length){
  if(feature == 'site'){
    take.col <- '$1,$2,$2'
  } else{
    take.col <- '$1,$2,$3'
  }
  x <- paste('awk -F"\\t" \'{print ', take.col, '}\' OFS="\\t" ',
             path.to.feature.stat,
             ' | tail -n +2 | bedtools slop -i stdin -g ',
             path.to.chrom.file,
             ' -b ', flank.length, ' | bedtools sort -i stdin| bedtools intersect -a ',
             path.to.ph.bed12,
             ' -b stdin -wa -wb', sep = '')
  df.handle <- read.delim(pipe(x), header = F)
  df.handle
}

IntersectPanhandleOne <- function(path.to.ph.bed12, path.to.feature.stat, feature){
  if(feature == 'site'){
    take.col <- '$1,$1,$2'
  } else{
    take.col <- '$1,$2,$3'
  }
  x <- paste('awk -F"\\t" \'{print $1,$2,$9,$4}\' OFS="\\t" ',
             path.to.ph.bed12, '> ph.bed')
  system(x)
  x <- paste('awk -F"\\t" \'{print ', take.col, '}\' OFS="\\t" ',
             path.to.feature.stat,
             ' | tail -n +2 ',
             ' | bedtools sort -i stdin| bedtools intersect -a ',
             'ph.bed',
             ' -b stdin -wa -wb', sep = '')
  df.ph <- read.delim(pipe(x), header = F)
  df.ph
}

IntersectPanhandles <- function(df, path.to.ph.main, path.to.feature.stat, flank.length, feature, loop.out = T, intron.retention = T){
  ph <- fread(path.to.ph.main)
  ph <- ph[ph$always.intronic == T & 
                   ph$CDS.transcript == T & 
                   ph$tRNA == F & 
                   ph$sno.miRNA == F & 
                   ph$TFBS == F & 
                   ph$Repeats == F, ]
  if(feature == 'site'){
    df.handle.left <- IntersectHandle(path.to.ph.bed12, path.to.chrom.file, path.to.feature.stat, feature, flank.length)
    SwapHandles(path.to.ph.bed12)
    df.handle.right <- IntersectHandle('./swapt_handles.bed', path.to.chrom.file, path.to.feature.stat, feature, flank.length)
    ph.left <- merge(ph, df.handle.left[, c('V14', 'V4')], by.x = 'id', by.y = 'V4', all.x = T)
    ph.right <- merge(ph, df.handle.right[, c('V14', 'V4')], by.x = 'id', by.y = 'V4', all.x = T)
    ph <- rbind(ph.left, ph.right)
    ph <- setnames(ph, old = c("V14"), new = c("feature_coord"))
    ph$feature_coord <- ph$feature_coord + flank.length
    ph <- merge(ph, df[, c('chr', 'coord', 'strand', 'passed.cutoff', 'delta_psi_c'), with = F], 
                by.x = c('chr', 'feature_coord', 'strand'),
                by.y = c('chr', 'coord', 'strand'), all.x = T)
  } else if(feature == 'intron'){
    if(intron.retention){
      df.ph <- IntersectPanhandleOne(path.to.ph.bed12, path.to.feature.stat, feature)
      df.ph <- df.ph[df.ph$V6 <= df.ph$V2 & df.ph$V3 <= df.ph$V7, ]
      df.ph <- df.ph[df.ph$V2 - df.ph$V6 <= flank.length & df.ph$V7 - df.ph$V3 <= flank.length, ]
      ph <- merge(ph, df.ph[, c('V4', 'V6', 'V7')], by.x = 'id', by.y = 'V4', all = F)
      ph <- setnames(ph, old = c("V6", "V7"), new = c("feature_start", "feature_end"))
    } else{
      df.handle.left <- IntersectHandle(path.to.ph.bed12, path.to.chrom.file, path.to.feature.stat, feature, flank.length)
      SwapHandles(path.to.ph.bed12)
      df.handle.right <- IntersectHandle('./swapt_handles.bed', path.to.chrom.file, path.to.feature.stat, feature, flank.length)
      ph.left <- merge(ph, df.handle.left[, c('V14', 'V15', 'V4')], by.x = 'id', by.y = 'V4')
      ph.right <- merge(ph, df.handle.right[, c('V14', 'V15', 'V4')], by.x = 'id', by.y = 'V4')
      ph <- rbind(ph.left, ph.right)
      ph <- setnames(ph, old = c("V14", "V15"), new = c("feature_start", "feature_end"))
      ph$feature_start <- ph$feature_start + flank.length
      ph$feature_end <- ph$feature_end - flank.length    
    }
  
    ph <- merge(ph, df[, c('chr', 'start', 'end', 'strand', 'passed.cutoff', 'delta_psi_c'), with = F], 
                by.x = c('chr', 'feature_start', 'feature_end', 'strand'),
                by.y = c('chr', 'start', 'end', 'strand'), all.x = T)

  } else if(feature == 'exon') {
    if(loop.out == T){
      df.ph <- IntersectPanhandleOne(path.to.ph.bed12, path.to.feature.stat, feature)
      df.ph <- df.ph[df.ph$V2 <= df.ph$V6 & df.ph$V7 <= df.ph$V3, ]
      ph <- merge(ph, df.ph[, c('V4', 'V6', 'V7')], by.x = 'id', by.y = 'V4', all = F)
      ph <- setnames(ph, old = c("V6", "V7"), new = c("feature_start", "feature_end"))
    } else{
      df.handle.left <- IntersectHandle(path.to.ph.bed12, path.to.chrom.file, path.to.feature.stat, feature, flank.length)
      SwapHandles(path.to.ph.bed12)
      df.handle.right <- IntersectHandle('./swapt_handles.bed', path.to.chrom.file, path.to.feature.stat, feature, flank.length)
      ph.left <- merge(ph, df.handle.left[, c('V14', 'V15', 'V4')], by.x = 'id', by.y = 'V4')
      ph.right <- merge(ph, df.handle.right[, c('V14', 'V15', 'V4')], by.x = 'id', by.y = 'V4')
      ph <- rbind(ph.left, ph.right)
      ph <- setnames(ph, old = c("V14", "V15"), new = c("feature_start", "feature_end"))
      ph$feature_start <- ph$feature_start + flank.length
      ph$feature_end <- ph$feature_end - flank.length
    }
    ph <- merge(ph, df[, c('chr', 'start', 'end', 'strand', 'passed.cutoff', 'delta_psi_c'), with = F], 
                by.x = c('chr', 'feature_start', 'feature_end', 'strand'),
                by.y = c('chr', 'start', 'end', 'strand'), all.x = T)
  }
  ph
}

MakeHandlesBed <- function(dt){
  dt.left <- dt[, c('chr', 'panhandle_start', 'panhandle_left_hand', 'id')]
  dt.right <- dt[, c('chr', 'panhandle_right_hand', 'panhandle_end', 'id')]
  colnames(dt.left) <- c('chr', 'start', 'stop', 'name')
  colnames(dt.right) <- c('chr', 'start', 'stop', 'name')
  dt.left$name <- paste(dt.left$name, '_left', sep = '')
  dt.right$name <- paste(dt.right$name, '_right', sep = '')
  dt.both <- rbind(dt.left, dt.right)
  dt.both
}

CountMutations <- function(handles, path.to.mut, title){
  x <- paste0('tail -n +2 ', handles, ' | bedtools sort -i stdin | bedtools merge | awk -F"\\t" \'{print $1,$3-$2+1}\' OFS="\\t"|awk \'{s+=$2}END{print s}\'')
  length_ <- tryCatch(scan(pipe(x)), error = function(e) 0)
  length_ <- as.numeric(length_)
  print(length_)
  x <- paste0('tail -n +2 ', handles, ' | bedtools sort -i stdin | bedtools merge | bedtools intersect -a stdin -b ', path.to.mut, ' -wa -wb')
  mut.in.df <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  if(is.null(mut.in.df)){
    count.mut <- 0
  } else{
    mut.in.df$mut.id <- str_split_fixed(mut.in.df$V7, '_', Inf)[, 1]
    if(!grepl('not_sign', handles)){
      write.table(mut.in.df, paste0(title, '_with_mut.tsv'), sep = '\t', row.names = F, quote = F)
    }
    count.mut <- length(unique(mut.in.df$mut.id))  
  }
  print(count.mut)
  return(c(length_, count.mut))
}

TestMutEnrichment <- function(df.ph, path.to.mut, title){
  ph.handles.sign <- MakeHandlesBed(df.ph[df.ph$passed.cutoff == T, ])
  ph.handles.not.sign <- MakeHandlesBed(df.ph[df.ph$passed.cutoff == F, ])
  write.table(ph.handles.sign, './handles_sign.bed', sep = '\t', row.names = F, quote = F)
  write.table(ph.handles.not.sign, './handles_not_sign.bed', sep = '\t', row.names = F, quote = F)
  
  out <- CountMutations('./handles_sign.bed', path.to.mut, title)
  length.sign <- out[1]
  count.mut.in.sign <- out[2]
  
  out <- CountMutations('./handles_not_sign.bed', path.to.mut, title)
  length.not.sign <- out[1]
  count.mut.in.not.sign <- out[2]
  
  p.value <- poisson.test(x = c(count.mut.in.sign, count.mut.in.not.sign), T = c(length.sign, length.not.sign), alternative = 't', r = 1)$p.value
  return(list('length.sign' = length.sign, 
              'count.mut.in.sign' = count.mut.in.sign,
              'length.not.sign' = length.not.sign,
              'count.mut.in.not.sign' = count.mut.in.not.sign,
              'p.value' = p.value))
}

IntersectWithMutatoins <- function(df.ph, path.to.mut, title){
  ph.handles.sign <- MakeHandlesBed(df.ph[df.ph$passed.cutoff == T, ])
  write.table(ph.handles.sign, './handles_sign.bed', sep = '\t', row.names = F, quote = F)
  x <- paste0('tail -n +2 ', handles, ' | bedtools sort -i stdin | bedtools intersect -a stdin -b ', path.to.mut, ' -wa -wb')
  mut.in.df <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  mut.in.df$id <- as.numeric(str_split_fixed(mut.in.df$V4, '_', Inf)[, 1])
  df.ph <- df.ph[df.ph$passed.cutoff == T, ]
  merged <- merge(df.ph, mut.in.df[, c('id', 'V4', 'V6', 'V7', 'V8')], by = 'id', all = F)
  setnames(merged, old = c('V4', 'V6', 'V7', 'V8'), new = c('handle.w.mut', 'mut_start', 'mut_end', 'mut_id'))
  
}

SelectCodingExons <- function(df, path.to.CDS.exons){
  write.table(df, 'tmp.bed', sep = '\t', col.names = F, row.names = F, quote = F)
  columns <- colnames(df)
  x <- paste0('bedtools sort -i tmp.bed | bedtools intersect -a stdin -b ', path.to.CDS.exons, ' -wa -wb')
  df <- as.data.table(read.delim(pipe(x), header = F))
  df <- df[, 1:(ncol(df) - 4)]
  colnames(df) <- columns
  df <- df[! duplicated(df), ]
  x <- 'rm tmp.bed'
  system(x)
  df
}

SaveTrack <- function(title, df.ph, feature){
  mut <- fread(paste0(title, '_with_mut.tsv'))
  mut <- mut[! duplicated(mut[, c('V1', 'V2', 'V3'), with = F]), ]
  to.save <- merge(df.ph, mut[, c('V1', 'V2', 'V3', 'V5', 'V7'), with = F], by.x = 'chr', by.y = c('V1'), all = F, allow.cartesian = TRUE)
  to.save <- to.save[to.save$panhandle_start <= to.save$V2 & to.save$V3 <= to.save$panhandle_left_hand | 
                       to.save$panhandle_left_hand <= to.save$V2 & to.save$V3 <= to.save$panhandle_end, ]
  to.save <- to.save[to.save$passed.cutoff == T, ]
  to.save <- to.save[! duplicated(to.save), ]
  to.save <- setnames(to.save, old = c('V5', 'V7'), new = c('mut_coord', 'mut_id'))
  if(feature == 'site'){
    to.save$index <- paste0(to.save$id, '_', to.save$energy, '_', to.save$feature_coord, '_', to.save$delta_psi_c)  
  } else{
    to.save$index <- paste0(to.save$id, '_', to.save$energy, '_', to.save$feature_start, '_', to.save$feature_end, '_', to.save$delta_psi_c)  
  }
  Make_bed(as.data.frame(to.save), paste0(title, '_track.bed'), different_colours = T, title)
}

main <- function(input, feature, stat, path.to.mut, title){
  main.results <- list()
  df <- fread(input)
  df <- Preprocess(df)
  if(feature == 'exon'){
    df <- SelectCodingExons(df, path.to.CDS.exons)
  }
  df <- CalculateSJcount(df, feature)
  df <- CalculatePsi(df)
  df <- df[df$sjcount_normal >= denom.cutoff & df$sjcount_tumor >= denom.cutoff, ]
  df <- RegressOutExpression(df)
  if(stat == 'qval'){
    #print('Calculating qval..')
    df <- IntersectWithGenes(df, denom.cutoff)
    df.stat <- CalculateQval(df, q.cutoff, p.cutoff, denom.cutoff)
    df$passed.cutoff <- F
    df.stat.passed <- df.stat[df.stat$q < q.cutoff, ]
    if (feature == 'site'){
      #print('use site..')
      df.stat.passed <- df.stat.passed[, c('chr', 'coord', 'strand', 'q'), with = F]
      df <- merge(df, df.stat.passed, by = c('chr', 'coord', 'strand'), all.x = T)  
    } else{
      #print('use intron or exon..')
      df.stat.passed <- df.stat.passed[, c('chr', 'start', 'end', 'strand', 'q'), with = F]
      df <- merge(df, df.stat.passed, by = c('chr', 'start', 'end', 'strand'), all.x = T)
    }
    df[!is.na(df$q), ]$passed.cutoff <- T 
  } else if(stat == 'delta.psi'){
    #print('Calculating delta psi..')
    df$passed.cutoff <- abs(df$delta_psi_c) > delta.psi.cutoff
  } else{
    print('Do not know such stat!')
  }
  print('Features')
  if(feature == 'site'){
    r1 <- dim(df[! duplicated(df[, c('chr' ,'coord','strand'), with = F]) & df$passed.cutoff, ])[1]
    print(dim(df[! duplicated(df[, c('chr' ,'coord','strand'), with = F]) & df$passed.cutoff, ])[1])
    
  } else{
    r1 <- dim(df[! duplicated(df[, c('chr' ,'start', 'end','strand'), with = F]) & df$passed.cutoff, ])[1]
    print(dim(df[! duplicated(df[, c('chr' ,'start', 'end','strand'), with = F]) & df$passed.cutoff, ])[1])
    
  }
  write.table(df, gsub('.tsv', '_stat.tsv', input), sep = '\t', row.names = F, quote = F)
  df.ph <- IntersectPanhandles(df, path.to.ph.main, 
                               path.to.feature.stat = gsub('.tsv', '_stat.tsv', input), 
                               flank.length, feature, loop.out = loop.out)
  print('Features with panhandles')
  if(feature == 'site'){
    r2 <- dim(df.ph[! duplicated(df.ph[, c('chr' ,'feature_coord','strand'), with = F]) & df.ph$passed.cutoff, ])[1]
    print(dim(df.ph[! duplicated(df.ph[, c('chr' ,'feature_coord','strand'), with = F]) & df.ph$passed.cutoff, ])[1])
    
  } else{
    r2 <- dim(df.ph[! duplicated(df.ph[, c('chr' ,'feature_start', 'feature_end','strand'), with = F]) & df.ph$passed.cutoff, ])[1]
    print(dim(df.ph[! duplicated(df.ph[, c('chr' ,'feature_start', 'feature_end','strand'), with = F]) & df.ph$passed.cutoff, ])[1])
    
  }
  results <- TestMutEnrichment(df.ph, path.to.mut, title)
  print(paste0('pvalue = ', as.character(round(results$p.value, 3)), 
               ' signif = ', as.character(round(results$count.mut.in.sign / results$length.sign, 5)), 
               ' not signif = ', as.character(round(results$count.mut.in.not.sign / results$length.not.sign, 5))))  
  main.results[['r1']] <- r1
  main.results[['r2']] <- r2
  main.results[['r3']] <- c(round(results$p.value, 3), 
                            round(results$count.mut.in.sign / results$length.sign, 5), 
                            round(results$count.mut.in.not.sign / results$length.not.sign, 5))
  
  SaveTrack(title, df.ph, feature)
  return(main.results)
}

ParseResults <- function(results){
  results2 <- list()
  for(exp in names(results)){
    results2[[exp]] <- c(results[[exp]]$r1, results[[exp]]$r2, results[[exp]]$r3)
  }
  results.df <- as.data.frame(t(as.data.frame(results2)))
  colnames(results.df) <- c('n.sign.features', 'n.sign.features.w.ph', 'pval', 'proc.sign', 'proc.not.sign')
  results.df$experiment <- rownames(results.df)
  results.df$cutoff.by <- str_split_fixed(results.df$experiment, '_', Inf)[, 1]
  results.df$feature <- str_split_fixed(results.df$experiment, '_', Inf)[, 2]
  results.df$tissue <- str_split_fixed(results.df$experiment, '_', Inf)[, 3]
  results.df <- results.df[, ! names(results.df) %in% c('experiment')]
  write.table(results.df, paste0(q.cutoff, '_', flank.length, '_results.tsv'), sep = '\t', row.names = F, quote = F)
}


results <- list()
for (tissue in c('kidney', 'liver')){
  for (stat in c('qval', 'delta.psi')){
    for (feature in c('site', 'exon', 'intron')){
      title <- paste0(stat, '_', feature, '_', tissue)
      print(title)
      input <- paste0('./data/', tissue, '_', feature, '.tsv')
      if(tissue == 'kidney'){
        path.to.mut <- '../mutations_icgc/SNVs/RECA_simple_somatic_mutation.open_only_SNV.bed'
      } else{
        path.to.mut <- '../mutations_icgc/SNVs/LIRI_simple_somatic_mutation.open_only_SNV_hypermut_donot_removed.bed'
      }
      results[[title]] <- main(input, feature, stat, path.to.mut, title)
    }
  }  
}
ParseResults(results)





p = ggplot(dt, aes(x = log10sjcount, y = delta_psi_c)) +
  geom_point(size = 0.5, alpha = 0.5, aes(color = (q < q.cutoff))) + 
  geom_abline(slope = 0,lty = "dashed") + 
  xlab("log(denomTumor + denomNormal)") +
  ylab(expression("dPSI(Tumor-Normal)")) + 
  theme_classic() + theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

png('~/Desktop/1.png', width = 6, height = 7, units = 'in', res = 300)
print(p)
dev.off()
