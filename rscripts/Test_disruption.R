#!/usr/bin/Rscript
library(stringr)
library(data.table)
library(ggplot2)
library(reshape2)
options(scipen = 999)

path.to.ph.w.mut <- "../python_scripts/folding_pretty_copy/out/out_mutator_kidney.tsv"
path.to.ph <- "../python_scripts/folding_pretty_copy/out/panhandles_preprocessed.tsv"
title <- 'kidney'

# Test_disruption.R path.to.ph.w.mut path.to.ph title 
args <- commandArgs(trailingOnly = TRUE)
path.to.ph.w.mut <- args[1]
path.to.ph <- args[2]
title <- args[3]

AnalyzeDeltaDeltaG <- function(dt, title){
  dt$dE.real <- dt$new_energy - dt$old_energy 
  dt$dE.random <- dt$random_energy - dt$old_energy 
  dt$delta.dE <- dt$dE.real - dt$dE.random
  png(paste0('./pictures/', title, '_delta_dE_boxplot.png'), 
      res = 300, width = 4, height = 5, units = 'in')
  medians <- data.table(median = median(dt$delta.dE))
  pvalue <- formatC(wilcox.test(dt$delta.dE, mu = 0, correct = F)$p.value, format = "e", digits = 1)
  print(ggplot(dt, aes(y = delta.dE)) + 
    geom_boxplot() +
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

dt <- fread(path.to.ph.w.mut)
dt <- AnalyzeDeltaDeltaG(dt, title)
#dt.melt <- AnalyzeDeltaG(dt, title)
#ph <- fread(path.to.ph)
# MakeTrack(dt, ph)








