options(scipen = 999)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gplots)
library(VennDiagram)
library(gridExtra)
library(ggpubr)

features <- c('exon', 'intron', 'site')
tissues <- c('kidney', 'liver')


pdf('data/compare_ipsa_spladder_scatterplot.pdf')
for(tissue in tissues){
  feature <- 'exon'
  spladder <- fread(paste0('data/spladder_output/spladder_', tissue, '_', feature, '_filtered_population.tsv'))
  ipsa <- fread(paste0('data/ipsa_output/ipsa_', tissue, '_', feature, '_filtered_population.tsv'))
  spladder$id <- paste(spladder$chr, spladder$start, spladder$end, sep = '_')
  ipsa$id <- paste(ipsa$chr, ipsa$start, ipsa$end, sep = '_')
  merged <- merge(spladder[, c('id', 'psi_normal'), with = F], ipsa[, c('id', 'psi_normal'), with = F], by = 'id', all = F, suffixes = c('_spladder', '_ipsa'))
  
  sp <- ggscatter(merged, x = "psi_normal_spladder", y = "psi_normal_ipsa",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  xlim = c(0, 1),
                  ylim = c(0, 1.1))
  # Add correlation coefficient
  sp <- sp + stat_cor(method = "pearson", label.x = 0, label.y = 1.1)
  print(sp)
}
dev.off()


pdf('data/compare_ipsa_spladder.pdf')
for(tissue in tissues){
  for(feature in features){
    print(paste0(tissue, '_', feature))
    if(feature == 'site'){
      spladder1 <- fread(paste0('data/spladder_output/spladder_', tissue, '_', '3prime_site', '_filtered_population.tsv'))
      spladder2 <- fread(paste0('data/spladder_output/spladder_', tissue, '_', '5prime_site', '_filtered_population.tsv'))
      spladder <- rbind(spladder1, spladder2)
      ipsa <- fread(paste0('data/ipsa_output/ipsa_', tissue, '_', feature, '_filtered_population.tsv'))
      spladder$id <- paste(spladder$chr, spladder$start, spladder$end, sep = '_')
      spladder$start_id <- paste(spladder$chr, spladder$start, sep = '_')
      spladder$end_id <- paste(spladder$chr, spladder$end, sep = '_')
      ipsa$id <- paste(ipsa$chr, ipsa$start, sep = '_')
      area1 <- length(unique(spladder$id))
      area2 <- length(unique(ipsa$id))
      cross.area <- length(unique(c(intersect(spladder$start_id, ipsa$id), intersect(spladder$end_id, ipsa$id))))
      
      area1.sign <- length(unique(spladder[spladder$significant, ]$id))
      area2.sign <- length(unique(ipsa[ipsa$significant, ]$id))
      cross.area.sign <- length(unique(c(intersect(spladder[spladder$significant, ]$start_id, ipsa[ipsa$significant, ]$id), 
                                         intersect(spladder[spladder$significant, ]$end_id, ipsa[ipsa$significant, ]$id))))
    } else{
      spladder <- fread(paste0('data/spladder_output/spladder_', tissue, '_', feature, '_filtered_population.tsv'))
      ipsa <- fread(paste0('data/ipsa_output/ipsa_', tissue, '_', feature, '_filtered_population.tsv'))
      spladder$id <- paste(spladder$chr, spladder$start, spladder$end, sep = '_')
      ipsa$id <- paste(ipsa$chr, ipsa$start, ipsa$end, sep = '_')
      area1 <- length(unique(spladder$id))
      area2 <- length(unique(ipsa$id))
      cross.area <- length(unique(intersect(spladder$id, ipsa$id)))
      
      area1.sign <- length(unique(spladder[spladder$significant, ]$id))
      area2.sign <- length(unique(ipsa[ipsa$significant, ]$id))
      cross.area.sign <- length(unique(intersect(spladder[spladder$significant, ]$id, ipsa[ipsa$significant, ]$id)))
    }
    print(paste(area1, area2, cross.area, area1.sign, area2.sign, cross.area.sign, sep = ', '))
    g <- draw.pairwise.venn(area1 = area1, 
                            area2 = area2, 
                            cross.area = cross.area, 
                            category = c("spladder", "ipsa"), 
                            fill = c("blue", "red"), 
                            alpha = rep(0.2, 2), 
                            cat.pos = c(0, 0), cat.dist = rep(0.025, 2),
                            lty = rep("blank", 2), ind = FALSE)
    grid.arrange(gTree(children = g), top = paste(tissue, feature, 'all', sep = ', '))
    
    g <- draw.pairwise.venn(area1 = area1.sign, 
                            area2 = area2.sign, 
                            cross.area = cross.area.sign, 
                            category = c("spladder", "ipsa"), 
                            fill = c("blue", "red"), 
                            alpha = rep(0.5, 2), 
                            cat.pos = c(0, 0), cat.dist = rep(0.025, 2),
                            lty = rep("blank", 2), ind = FALSE)
    grid.arrange(gTree(children = g), top = paste(tissue, feature, 'significant', sep = ', '))    
  }
}
dev.off()

