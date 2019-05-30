options(scipen = 999)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(grid)
library(gtable)

MakeFlankedHandles <- function(ph.handle.path.bed, flank.length = 1000){
  handles <- fread(ph.handle.path.bed)
  handles$V2 <- handles$V2 - flank.length
  handles$V3 <- handles$V3 + flank.length
  write.table(handles, gsub('.bed', '_with_flanks.bed', ph.handle.path.bed), row.names = F, col.names = F, quote = F, sep = '\t')
}

IntersectPanhandles <- function(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out){
  feature.path <- paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')
  print(feature.path)
  feature.dt <- fread(feature.path)
  if(feature == 'exon' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.path.bed, ' -wa > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)    
  } else if(feature == 'exon' & looped.out == F){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)   
  } else if(grepl('site', feature)){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)   
  }
  else if(feature == 'intron' & looped.out == F){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x) 
    dt <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
    ncol <- dim(feature.dt)[2]
    dt <- dt[dt$V2 <= dt[[ncol + 3]] & dt$V2 >= dt[[ncol + 2]] | dt$V3 <= dt[[ncol + 3]] & dt$V3 >= dt[[ncol + 2]], ]
    write.table(dt, paste0(prefix, tissue, '_', feature, '_overlapped.tsv'), row.names = F, col.names = F, quote = F, sep = '\t')
  } else if(feature == 'intron' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x) 
    dt <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
    ncol <- dim(feature.dt)[2]
    dt <- dt[dt$V2 <= dt[[ncol + 3]] & dt$V2 >= dt[[ncol + 2]] | dt$V3 <= dt[[ncol + 3]] & dt$V3 >= dt[[ncol + 2]], ]
    dt$id <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 1]
    dt$handle <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 2]
    overlapped <- merge(dt, dt, by = c(paste("V", c(1:ncol), sep = ''), 'id'))
    overlapped <- overlapped[overlapped$handle.y != overlapped$handle.x, ]
    overlapped$V21.x <- overlapped$V21.x + 1000
    overlapped$V21.y <- overlapped$V21.y + 1000
    overlapped$V22.x <- overlapped$V22.x - 1000
    overlapped$V22.y <- overlapped$V22.y - 1000
    overlapped <- overlapped[overlapped$V2 <= overlapped$V21.x & 
                               overlapped$V2 <= overlapped$V21.y & 
                               overlapped$V3 >= overlapped$V22.x & 
                               overlapped$V3 >= overlapped$V22.y,]
    to.drop <- grep('x|y', colnames(overlapped), value = T)
    overlapped <- overlapped[, ! names(overlapped) %in% c('id', to.drop), with = F]
    overlapped <- overlapped[! duplicated(overlapped[, c('V1', "V2", 'V3', "V4"), with = F]), ]
    write.table(overlapped, paste0(prefix, tissue, '_', feature, '_overlapped.tsv'), row.names = F, col.names = F, quote = F, sep = '\t')
  }
  overlapped <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
  overlapped$overlap <- T
  dt <- merge(feature.dt, overlapped[, c('V1', 'V2', 'V3', 'V4', 'overlap'), with = F], 
              by.x = c('chr', 'start', 'end', 'strand'), 
              by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  dt[is.na(dt$overlap), ]$overlap <- F
  dt <- dt[!duplicated(dt[,c('chr','start','end', 'strand'), with = F]), ]
  if(grepl('site', feature)){
    dt <- dt[dt$end - dt$start <= 500, ]
  }
  dt
}

PlotChisq <- function(df, title){
  title <- textGrob(title, gp = gpar(fontsize = 10))
  padding <- unit(0.5, "line")
  df <- tableGrob(df)
  table <- gtable_add_rows(
    df, heights = grobHeight(title) + padding, pos = 0
  )
  table <- gtable_add_grob(
    table, list(title),
    t = 1, l = 1, r = ncol(table)
  )
  grid.newpage()
  grid.draw(table)
}

TestPsi <- function(tissue, feature, locations, ph.path.bed, ph.handle.flanked.path.bed, out.path, prefix, exons.bridges.path){
  counts <- list()
  dt.nlooped.out <- IntersectPanhandles(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out = F)
  dt.nlooped.out <- dt.nlooped.out[, c('chr', 'start', 'end', 'strand', 'psi_normal', 'psi_tumor', 'delta_psi', 'delta_psi_c', 'significant', 'overlap'), with = F]
  dt.nlooped.out$feature.id <- paste(dt.nlooped.out$chr, dt.nlooped.out$start, dt.nlooped.out$end, sep = '_')
  if(feature == 'exon' & looped.out == F){
    exons.bridges <- fread(exons.bridges.path)
    exons.bridges$feature_id <- paste(exons.bridges$chr, exons.bridges$exon_start, exons.bridges$exon_end, sep = '_')
    dt.nlooped.out$eCLIP <- F
    dt.nlooped.out[dt.nlooped.out$feature.id %in% exons.bridges$feature_id, ]$eCLIP <- T
  }
  counts[['total_denom_filtered']] <- dim(dt.nlooped.out)[1]
  counts[['DE']] <- dim(dt.nlooped.out[dt.nlooped.out$significant, ])[1]
  counts[['close_struct']] <- dim(dt.nlooped.out[dt.nlooped.out$overlap, ])[1]
  counts[['close_struct_and_DE']] <- dim(dt.nlooped.out[dt.nlooped.out$overlap & dt.nlooped.out$significant, ])[1]
  
  if(length(locations) > 1){
    dt.looped.out <- IntersectPanhandles(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out = T)
    dt.looped.out <- dt.looped.out[, c('chr', 'start', 'end', 'strand', 'psi_normal', 'psi_tumor', 'delta_psi', 'delta_psi_c', 'significant', 'overlap'), with = F]
    counts[['looped.out']] <- dim(dt.looped.out[dt.looped.out$overlap, ])[1]
    counts[['looped.out_and_DE']] <- dim(dt.looped.out[dt.looped.out$overlap & dt.looped.out$significant, ])[1]
  }
  colours <- c('#F8766D', '#00BFC4')
  names(colours) <- c(TRUE, FALSE)
  pdf(paste0(out.path, '/', tissue, '_', feature, '.pdf'))
  for(looped.out in locations){
    ### PSI
    if(looped.out){
      dt <- dt.looped.out
      if(dim(dt[dt$overlap == T, ])[1] > 0 & dim(dt[dt$overlap == F, ])[1] > 0 ){
        p.value <- wilcox.test(dt[dt$overlap == T, ]$psi_normal, dt[dt$overlap == F, ]$psi_normal, correct = F)$p.value
      } else {
        p.value <- 'NA'
      }
      g <- ggplot(dt, aes(x = psi_normal, fill = overlap)) +
        geom_histogram(position = 'dodge', aes(y=..density..)) +
        xlab('PSI') +
        ggtitle(label = paste0(tissue, ', ', feature), subtitle = paste0('p-value = ', formatC(p.value, format = "e", digits = 2))) +
        scale_fill_manual(values = colours) +
        guides(fill = guide_legend(title = "looped out"))
      print(g)
      means <- aggregate(psi_normal ~  overlap, dt, mean)
      means$psi_normal <- round(means$psi_normal, 2)
      g <- ggplot(dt, aes(x = overlap, y = psi_normal, fill = overlap)) +
        geom_boxplot() + 
        scale_fill_manual(values = colours) +
        ggtitle(label = paste0(tissue, ', ', feature), subtitle = paste0('p-value = ', formatC(p.value, format = "e", digits = 2))) +
        guides(fill = guide_legend(title = "looped out")) + 
        stat_summary(fun.y = mean, colour = "darkred", geom = "point", shape = 18, size = 3, show_guide = FALSE) + 
        geom_text(data = means, aes(label = psi_normal, y = psi_normal + 0.02))
        
      print(g)
      ### delta PSI
      if(dim(dt[dt$significant & dt$overlap == T, ])[1] > 0 & dim(dt[dt$significant & dt$overlap == F, ])[1] > 0 ){
        # p.value <- wilcox.test(dt[dt$significant & dt$overlap == T, ]$delta_psi_c, 
        #                        dt[dt$significant & dt$overlap == F, ]$delta_psi_c, correct = F, alternative = 'l')$p.value
        # p.value2 <- wilcox.test(dt[dt$significant & dt$overlap == T, ]$delta_psi_c, 
        #                        dt[dt$significant & dt$overlap == F, ]$delta_psi_c, correct = F, alternative = 't')$p.value
        p.value <- wilcox.test(x = dt[dt$overlap == T, ]$delta_psi_c, alternative = 'g', paired = F, correct = F)$p.value
        p.value2 <- wilcox.test(x = dt[dt$overlap == F, ]$delta_psi_c, alternative = 'g', paired = F, correct = F)$p.value
      } else {
        p.value <- 'NA'
        p.value2 <- 'NA'
      }
      means <- aggregate(delta_psi_c ~  overlap, dt[dt$significant], mean)
      means$delta_psi_c <- round(means$delta_psi_c, 3)
      g <- ggplot(dt[dt$significant, ], aes(x = overlap, y = delta_psi_c, fill = overlap)) +
        #geom_violin() +
        geom_boxplot() +
        stat_summary(fun.y = mean, colour = "darkred", geom = "point", shape = 18, size = 3, show_guide = FALSE) + 
        geom_text(data = means, aes(label = delta_psi_c, y = delta_psi_c - 0.1)) +
        xlab('') +
        ylab('delta PSI') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ggtitle(label = paste0(tissue, ', ', feature), subtitle = paste0('p-value = ', formatC(p.value2, format = "e", digits = 2), 
                                                                         '\tp-value = ', formatC(p.value, format = "e", digits = 2))) +
        scale_fill_manual(values = colours) +
        guides(fill = guides(fill = guide_legend(title = "looped out"))) +
        geom_hline(yintercept = 0, col = 'red')
      print(g)
      ### abs delta PSI
      if(dim(dt[dt$significant & dt$overlap == T, ])[1] > 0 & dim(dt[dt$significant & dt$overlap == F, ])[1] > 0 ){
        p.value <- wilcox.test(abs(dt[dt$significant & dt$overlap == T, ]$delta_psi_c), abs(dt[dt$significant & dt$overlap == F, ]$delta_psi_c), correct = F)$p.value
      } else {
        p.value <- 'NA'
      }
      g <- ggplot(dt[dt$significant, ], aes(x = overlap, y = abs(delta_psi_c), fill = overlap)) +
        geom_boxplot() +
        #geom_violin() +
        xlab('') +
        ylab('|delta PSI|') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ggtitle(label = paste0(tissue, ', ', feature), subtitle = paste0('p-value = ', formatC(p.value, format = "e", digits = 2))) +
        scale_fill_manual(values = colours)
      if(looped.out){
        g <- g + guides(fill = guides(fill = guide_legend(title = "looped out"))) 
      } else{
        g <- g + guides(fill = guides(fill = guide_legend(title = "close RNA structure"))) 
      }
      print(g)
    } else{
      ### abs delta PSI
      dt <- dt.nlooped.out
      if(dim(dt[dt$significant & dt$overlap, ])[1] > 0 & dim(dt[dt$significant & dt$overlap == F, ])[1] > 0 ){
        p.value <- wilcox.test(abs(dt[dt$significant & dt$overlap, ]$delta_psi_c), abs(dt[dt$significant & dt$overlap == F, ]$delta_psi_c), correct = F)$p.value
      } else {
        p.value <- 'NA'
      }
      g <- ggplot(dt[dt$significant, ], aes(x = overlap, y = abs(delta_psi_c), fill = overlap)) +
        geom_boxplot() +
        #geom_violin() +
        xlab('') +
        ylab('|delta PSI|') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ggtitle(label = paste0(tissue, ', ', feature), subtitle = paste0('p-value = ', formatC(p.value, format = "e", digits = 2))) +
        scale_fill_manual(values = colours)
      if(looped.out){
        g <- g + guides(fill = guides(fill = guide_legend(title = "looped out"))) 
      } else{
        g <- g + guides(fill = guides(fill = guide_legend(title = "close RNA structure"))) 
      }
      print(g)      
    }
  
    ##### chisq tests
    if(looped.out){
      dt <- dt.looped.out
      ### delta PSI 
      df <- data.frame(overlapped = c(dim(dt[dt$overlap & dt$delta_psi_c > 0, ])[1], dim(dt[dt$overlap & dt$delta_psi_c < 0, ])[1]),
                       not.overlapped = c(dim(dt[dt$overlap == F & dt$delta_psi_c > 0, ])[1], dim(dt[dt$overlap == F & dt$delta_psi_c < 0, ])[1]))
      row.names(df) <- c('delta PSI > 0', 'delta PSI < 0')
      PlotChisq(df, title =  paste0('LOOPED OUT \n all delta PSI > 0:\n overlapped = ', 100*round(dim(dt[dt$overlap & dt$delta_psi_c > 0, ])[1] / (dim(dt[dt$overlap & dt$delta_psi_c != 0, ])[1]), 2), '%\n',
                                             'not overlapped = ', 100*round(dim(dt[dt$overlap == F & dt$delta_psi_c > 0, ])[1] / (dim(dt[dt$overlap == F & dt$delta_psi_c != 0, ])[1]), 2),'%\n',
                                             'chisq test p-value = ', round(chisq.test(df, correct = T)$p.value, 2)))
      ### delta psi significant
      df <- data.frame(overlapped = c(dim(dt[dt$significant & dt$overlap & dt$delta_psi_c > 0, ])[1], dim(dt[dt$significant & dt$overlap & dt$delta_psi_c < 0, ])[1]),
                       not.overlapped = c(dim(dt[dt$significant & dt$overlap == F & dt$delta_psi_c > 0, ])[1], dim(dt[dt$significant & dt$overlap == F & dt$delta_psi_c < 0, ])[1]))
      row.names(df) <- c('delta PSI > 0', 'delta PSI < 0')
      PlotChisq(df, title =  paste0('LOOPED OUT \n significant delta PSI > 0: \n overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap & dt$delta_psi_c > 0, ])[1] / 
                                                                                                                      (dim(dt[dt$significant & dt$overlap & dt$delta_psi_c != 0, ])[1]), 2), '%\n',
                                             'not overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap == F & dt$delta_psi_c > 0, ])[1] / 
                                                                                                                      (dim(dt[dt$significant & dt$overlap == F & dt$delta_psi_c != 0, ])[1]), 2),'%\n',
                                             'chisq test p-value = ', round(chisq.test(df, correct = T)$p.value, 2)))
      
      ### significant
      df <- data.frame(overlapped = c(dim(dt[dt$significant & dt$overlap, ])[1], dim(dt[!dt$significant & dt$overlap, ])[1]),
                       not.overlapped = c(dim(dt[dt$significant & dt$overlap == F , ])[1], dim(dt[!dt$significant & dt$overlap == F , ])[1]))
      row.names(df) <- c('significant', 'not significant')
      PlotChisq(df, title =  paste0('LOOPED OUT \nsignificant: \n overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap , ])[1] / (dim(dt[dt$overlap, ])[1]), 3), '%\n',
                                             'not overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap == F, ])[1] / (dim(dt[dt$overlap == F, ])[1]), 3),'%\n',
                                             'chisq test p-value = ', round(chisq.test(df, correct = T)$p.value, 2)))

    } else{
      dt <- dt.nlooped.out
      ### significant
      df <- data.frame(overlapped = c(dim(dt[dt$significant & dt$overlap, ])[1], dim(dt[!dt$significant & dt$overlap, ])[1]),
                       not.overlapped = c(dim(dt[dt$significant & dt$overlap == F , ])[1], dim(dt[!dt$significant & dt$overlap == F , ])[1]))
      row.names(df) <- c('significant', 'not significant')
      PlotChisq(df, title =  paste0('CLOSE RNA \n significant: \n overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap , ])[1] / (dim(dt[dt$overlap, ])[1]), 3), '%\n',
                      'not overlapped = ', 100*round(dim(dt[dt$significant & dt$overlap == F, ])[1] / (dim(dt[dt$overlap == F, ])[1]), 3),'%\n',
                      'chisq test p-value = ', round(chisq.test(df, correct = T)$p.value, 2)))   

      g <- ggplot(dt, aes(x = eCLIP, y = abs(delta_psi_c), fill = eCLIP)) +
        geom_boxplot() +
        #geom_violin() +
        xlab('') +
        ylab('|delta PSI|') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ggtitle(label = paste0(tissue, ', ', feature)) +
        scale_fill_manual(values = colours)
      
      df <- data.frame(bridge = c(length(unique(dt[dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[!dt$significant & dt$eCLIP, ]$feature.id))), 
                        no.bridge = c(length(unique(dt[dt$significant & !dt$eCLIP, ]$feature.id)),length(unique(dt[!dt$significant & !dt$eCLIP, ]$feature.id))))
      row.names(df) <- c('significant', 'not significant')
      PlotChisq(df, title =  paste0('CLOSE RNA \n significant: \n bridge = ', 100*round(unname(unlist(df['significant',] / (df['significant',] + df['not significant',])))[1], 3), '%\n',
                                    'no bridge = ', 100*round(unname(unlist(df['significant',] / (df['significant',] + df['not significant',])))[2], 3),'%\n',
                                    'fisher test p-value = ', round(fisher.test(df)$p.value, 2)))   
      
      
      df <- data.frame(bridge = c(length(unique(dt[dt$significant & dt$eCLIP, ]$feature.id)),length(unique(dt[!dt$significant & dt$eCLIP, ]$feature.id))), 
                       no.bridge = c(length(unique(dt[dt$significant & !dt$eCLIP & dt$overlap, ]$feature.id)),length(unique(dt[!dt$significant & !dt$eCLIP & dt$overlap, ]$feature.id))))
      row.names(df) <- c('significant', 'not significant')
      PlotChisq(df, title =  paste0('CLOSE RNA \n significant: \n bridge = ', 100*round(unname(unlist(df['significant',] / (df['significant',] + df['not significant',])))[1], 3), '%\n',
                                    'ph = ', 100*round(unname(unlist(df['significant',] / (df['significant',] + df['not significant',])))[2], 3),'%\n',
                                    'fisher test p-value = ', round(fisher.test(df)$p.value, 2)))
    }

  }

  dev.off()  
  counts
}

exons.bridges.path <- '../RBPs/eCLIP_and_KD.tsv'
programmes <- c('spladder', 'ipsa')
datasets <- c('whole', 'filtered')
tissues <- c('liver', 'kidney')
for(programme in programmes){
  if(programme == 'spladder'){
    features <- c('exon', 'intron', '5prime_site', '3prime_site')
    out.path <- 'data/spladder_output/'
    prefix <- 'data/spladder_output/spladder_'
  } else if(programme == 'ipsa'){
    features <- c('exon', 'intron', 'site')
    out.path <- 'data/ipsa_output/'
    prefix <- 'data/ipsa_output/ipsa_'
  }
  for(dataset in datasets)
  {
    if(dataset == 'filtered'){
      ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered_handles.bed'
      ph.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered.bed'
      out.path2 <- paste0(out.path, 'filtered_ds_pictures/')
    } else if(dataset == 'whole'){
      ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_handles.bed'
      ph.path.bed <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed.bed'
      out.path2 <- paste0(out.path, 'whole_ds_pictures/')
    }
    ph.handle.flanked.path.bed <- gsub('.bed', '_with_flanks.bed', ph.handle.path.bed)
    if (!file.exists(ph.handle.flanked.path.bed)){
      MakeFlankedHandles(ph.handle.path.bed)
    }
    counts.list <- list()
    for(tissue in tissues){
      for(feature in features){
        if(grepl('site', feature)){
          locations <- c(FALSE)
        } else{
          locations <- c(FALSE, TRUE)
        }
        print(paste0(programme, '_', dataset, '_', tissue, '_', feature))
        counts.list[[paste0(tissue, '_', feature)]] <- TestPsi(tissue, feature, locations, ph.path.bed, ph.handle.flanked.path.bed, out.path2, prefix, exons.bridges.path)
      }
    }
    counts.table <- as.data.table(counts.list)
    for(column in colnames(counts.table)){
      set(counts.table, j = column, value = unlist(counts.table[[column]]))
    }
    write.table(counts.table, paste0(out.path2, '/', 'counts_table.tsv'), sep = '\t', row.names = F)
  }}









