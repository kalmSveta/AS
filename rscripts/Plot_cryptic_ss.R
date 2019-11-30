MakeCategories2 <- function(dt, what = 'energy', categories, strict = F){
  dt.new <- copy(dt)
  dt.new[, paste0(what, '.bin') := categories[1]]
  for(cat in categories[2:length(categories)]){
    if(strict){
      dt.tmp <- dt[dt[[what]] < cat, ]
    } else{
      dt.tmp <- dt[dt[[what]] <= cat, ]  
    }
    dt.tmp[, paste0(what, '.bin') := cat]
    dt.new <- rbind(dt.new, dt.tmp)
  }
  #dt.new <- dt.new[!duplicated(dt.new), ]
  dt.new
}


PlotTypesCutoffs <- function(counts.all, what. = 'introns', bin = 'energy.bin', colors. = energy.colors, shuffle){
  tmp <- merge(subset(counts.all, experiment.type == 'real'), 
               subset(counts.all, experiment.type != 'real'), 
               by = c('intersection.type', bin, 'expression'))
  tmp$ratio <- tmp$counts.x / tmp$counts.y
  title. <-  if(shuffle == 'phs') 'Random shift' else 'Random gene' 
  g <- ggplot(tmp, aes(y = ratio, x = intersection.type)) + 
    geom_boxplot(aes(col = get(bin))) +
    theme_linedraw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 15)) +
    xlab('') +
    ylab('enrichment') + 
    geom_hline(yintercept = 1, col = 'darkred', linetype = "dashed") +
    ggtitle(title.)
  if(bin == 'energy.bin'){
    g <- g + 
      scale_color_manual(values = colors., name = '\u0394G, kcal/mol')
  } else if(bin == 'spread.bim'){
    g <- g + 
      scale_color_manual(values = colors., name = 'spread, nts')
  } else if(bin == 'E-value.bin'){
    g <- g + 
      scale_color_manual(values = colors., name = 'E-value')
  }
  g <- g + facet_wrap(~expression, nrow = 1)
  g
}

ProcessSS <- function(dt, phs, evals, counts.cutoff = 9){
  phs.dt <- fread(phs)
  phs.dt$id <- str_split_fixed(phs.dt$V4, '_', 2)[, 1]
  phs.dt$energy <- as.numeric(str_split_fixed(phs.dt$V4, '_', 2)[, 2])
  dt$id <- as.character(str_split_fixed(dt$ph.id, '_', 2)[, 1])
  dt <- merge(dt, phs.dt[, c('id', 'energy'), with = F], by = 'id', all.x = T)
  evals.dt <- fread(evals)
  evals.dt$id <- as.character(evals.dt$id)
  dt <- merge(dt, evals.dt, by = 'id', all.x = T)
  dt$expression <- 'Low expression'
  dt[dt$counts >= counts.cutoff, ]$expression <- 'High expression'
  dt
}



DoIntesectioncutoffs <- function(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, db = 'counts'){
  dt <- fread(paste0(path.to.ph, 'phs_', what, '_shuffle_', shuffle, '_db_', db, '_position.tsv'))
  dt <- ProcessSS(dt, phs, evals, 10)
  
  energy.cutoffs <- c(-15, -20, -25, -30)
  energy.colors <- brewer.pal(n = 12, "Paired")[c(4, 7, 8, 6)]
  names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')
  column = "energy"
  dt.energy <- MakeCategories2(dt = dt, what = column, categories = energy.cutoffs, strict = F)
  dt.energy$energy.bin <- paste('<=', dt.energy$energy.bin, sep = '')
  dt.energy$energy.bin <- factor(dt.energy$energy.bin)
  dt2 <- dt.energy[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'energy.bin', 'expression'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + energy.bin + expression, data = dt2, FUN = sum))
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_', N.shuffle, '_db_', db, '_cutoff_energy.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'energy.bin', colors. = energy.colors, shuffle)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_db_', db, '_energy_cutoffs', '.pdf'))
  print(p)
  dev.off()
  
  eval.cutoffs <- c(2, 1, 0.05)
  eval.colors <- brewer.pal(n = 9, "Set1")[c(9, 5, 2)]
  names(eval.colors) <- paste('<', eval.cutoffs, sep = '')
  names(eval.colors)[1] <- 'all'
  column = "E-value"
  dt.eval <- MakeCategories2(dt = dt, what = column, categories = eval.cutoffs, strict = T)
  dt.eval$`E-value.bin` <- paste('<', dt.eval$`E-value.bin`, sep = '')
  dt.eval[dt.eval$`E-value.bin` == '<2']$`E-value.bin` <- 'all'
  dt.eval$`E-value.bin` <- factor(dt.eval$`E-value.bin`)
  dt2 <- dt.eval[, c('ph.id', 'counts', 'intersection.type', 'experiment.type', 'E-value.bin', 'expression'), with = F]
  counts.all <- as.data.table(aggregate(counts ~ experiment.type + intersection.type + `E-value.bin` + expression, data = dt2, FUN = sum))
  write.table(counts.all, paste0('intersection_', what,  '_shuffle_', shuffle, '_', N.shuffle, '_db_', db, '_cutoff_evalue.tsv'), sep = '\t', row.names = F)
  p <- PlotTypesCutoffs(counts.all, what. = what, bin = 'E-value.bin', colors. = eval.colors, shuffle)
  cairo_pdf(paste0('intersection_', what, '_shuffle_', shuffle, '_',  as.character(N.shuffle), '_db_', db, '_evalue_cutoffs', '.pdf'))
  print(p)
  dev.off()
  
}


path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
phs <- paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')
N.shuffle <- 10
genes <- '../conservative_features/not_intersected_coding_genes.bed'
conins <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.bed'
evals <- paste0(path.to.ph, 'R_scape_estended_all_pretty.tsv')

db = 'counts'
what = 'cryptic_ss'
intervals <- '../cryptic_ss/K562_cryptic_ss.bed'
shuffle = 'phs'
print(what)
print(shuffle)
DoIntesection(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph, conins)
DoIntesectioncutoffs(phs, genes, intervals, N.shuffle, what, shuffle, path.to.ph)   