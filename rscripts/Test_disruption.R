#!/usr/bin/Rscript
library(stringr)
library(data.table)
library(ggplot2)
library(reshape2)
options(scipen = 999)

path.to.ph.w.mut <- "../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/out_mutator_population_shuffle_nts_pos_new.tsv"
path.to.ph <- "../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv"
title <- 'population_pos_mut_filtered'

# Test_disruption.R path.to.ph.w.mut path.to.ph title 
args <- commandArgs(trailingOnly = TRUE)
path.to.ph.w.mut <- args[1]
path.to.ph <- args[2]
title <- args[3]

AnalyzeDeltaDeltaG <- function(dt, title){
  dt$dE.real <- dt$new_energy - dt$old_energy 
  dt$dE.random <- dt$random_energy - dt$old_energy 
  dt$delta.dE <- dt$dE.real - dt$dE.random
  dt <- dt[dt$delta.dE != 0, ]
  dt.small <- copy(dt)
  dt.small <- subset(dt, abs(delta.dE) > 2)
  #dt.small <- dt[dt$`E-value` < 0.01, ]
  dt.small$x <- 1
  pdf(gsub('.tsv', '_delta_dE_boxplot.pdf', path.to.ph.w.mut))
  medians <- round(data.table(median = median(dt.small$delta.dE)), 2)
  pvalue <- formatC(wilcox.test(dt.small$delta.dE, mu = 0, correct = F)$p.value, format = "e", digits = 1)
  print(ggplot(dt.small, aes(y = delta.dE)) + 
    geom_boxplot() +
    #geom_violin() + 
    theme_classic() +
    xlab(label = '') +
    ylim(-10, 10) +
    ylab(label = 'dG real - dG random, kcal/mol') +
    theme(axis.title = element_text(size = 20), 
          axis.text.y = element_text(size = 15), 
          #axis.line.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    geom_hline(yintercept = 0, col = 'red') + 
    geom_text(data = medians, aes(label = paste('median =', median, sep = ' '), y = median - 0.5, x = 0), size = 5) + 
    ggtitle(subtitle = paste0('p.value = ', pvalue), label = title))
  dev.off()
  dt
}

AnalyzeDeltaG <- function(dt, title){
  dt.melt <- dt[, c('mutation_info', 'panhandle_id', "dE.real", "dE.random")]
  dt.melt <- melt(dt.melt, id.vars = c('mutation_info', 'panhandle_id'))
  print('Wilcox test for dEreal vs dErandom:')
  print(wilcox.test(x = dt$dE.real, y = dt$dE.random, paired = F, mu = 0, correct = F))
  mean(dt$delta.dE)
  median(dt$delta.dE)
  png(paste0('./pictures/', title, '_delta_dE_boxplot_2on1.png'), 
      res = 300, width = 7, height = 5, units = 'in')
  #medians <- aggregate(value ~  variable, dt.melt, median)
  print(ggplot(dt.melt, aes(x = variable, y = value)) + 
    geom_boxplot() +
    theme_classic() +
    xlab(label = '') +
    ylab(label = 'dG, kcal/mol') +
    scale_x_discrete(labels = c("real", "random")) +
    theme(axis.title = element_text(size = 20), 
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 15)) +
    geom_hline(yintercept = 0, col = 'red')  + 
    ggtitle(label = title))
  #geom_text(data = medians, aes(label = paste('median =', value, sep = ' '), y = value - 0.5), size = 5)
  dev.off()
  dt.melt
}

MakeTrack <- function(dt, ph){
  merged <- merge(dt, ph, by.x = 'panhandle_id', by.y = 'id', all.x = T)
  merged$dE.real <- round(merged$dE.real, 2) 
  file_name <- gsub(".tsv", "_with_ph.tsv", path.to.ph.w.mut)
  merged$itemRGB <- '0,0,0'
  merged[merged$new_energy > -15, ]$itemRGB <- '255,0,0'
  write.table(merged, file_name, sep = '\t', row.names = F, quote = F)
  system(paste("./Make_big_bed.R", file_name, paste0("panahndles_with_mutations_from_", title), "panhandle_id,new_energy"))
}

dt <- fread(path.to.ph.w.mut, fill = T)
dt <- dt[dt$random_structure != '',]
dt <- AnalyzeDeltaDeltaG(dt, title)
#dt.melt <- AnalyzeDeltaG(dt, title)
#ph <- fread(path.to.ph)
# MakeTrack(dt, ph)



# qval <- fread('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/R_scape_estended_all_pretty.tsv')
# dt <- merge(dt, qval, by.x = 'panhandle_id', by.y = 'id')

ph <- fread(path.to.ph)
dt <- merge(dt, ph[, c('id', 'chr'), with = F], by.x = 'panhandle_id', by.y = 'id')

path.to.mut <- '../1000genomes/new/hg19_ss_flanks/phs/'
list.of.files <- file.info(list.files(path.to.mut, '*counts.2', full.names = T))
list.of.non.empty.files <- subset(list.of.files, size != 0)
mut.files <- rownames(list.of.non.empty.files)
dt.list <- lapply(mut.files, function(file){
  print(file)
  mut <- fread(file)
  mut <- subset(mut, V6 > 25)
  mut$V1 <- paste('chr', mut$V1, sep = '')
  dt2 <- merge(dt, mut[, c('V1','V2','V5','V6'), with = F], 
               by.x = c('chr', 'mut_coord', 'mut_to'), 
               by.y = c('V1', 'V2', 'V5'), all = F)
  subset(dt2, !is.na(dt2$V6))
})

dt <- Reduce(rbind, dt.list)
