.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(Cairo)
library(grid)
options(scipen = 999)

# Rscript FDR_plots.R path.to.ph.real path.to.ph.random folder.out

args <- commandArgs(trailingOnly = TRUE)
path.to.ph.real <- args[1]
path.to.ph.random <- args[2]
folder.out <- args[3]

# path.to.ph.real <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles.tsv'
# path.to.ph.random <- "../python_scripts/folding_pretty_copy/out/hg19_strand_specific/FDR/with_GC_content/"
# folder.out <- "../python_scripts/folding_pretty_copy/out/hg19_strand_specific/FDR/with_GC_content/"

MakeTable <- function(path.to.ph.random, path.to.ph.real, path.out, Xs, column = "energy"){
  dts <- list()
  for(file in list.files(path.to.ph.random, pattern = 'random_panhandles[0-9]*[0-9].tsv', full.names = T)){
    reals <- c()
    randoms <- c()
    ph.real <- fread(path.to.ph.real)
    ph.random <- fread(file)

    ph.real$interval1_start <- as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 2])
    ph.real$interval2_start <- as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 2])
    ph.real$loop.length <- ph.real$interval2_start + ph.real$start_al2 - (ph.real$interval1_start + ph.real$end_al1)
    
    ph.random$interval1_start <- as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 2])
    ph.random$interval2_start <- as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 2])
    ph.random$loop.length <- ph.random$interval2_start + ph.random$start_al2 - (ph.random$interval1_start + ph.random$end_al1)
    
    
    ph.real <- ph.real[ph.real$loop.length <= 10000, ]
    ph.random <- ph.random[ph.random$loop.length <= 10000, ]
    
    for(x in Xs){
      randoms <- c(randoms, dim(subset(ph.random, get(column) <= x))[1])
      reals <- c(reals, dim(subset(ph.real, get(column) <= x))[1])
    }
    dt <- data.table('cutoff' = Xs, 'FDR' = randoms / reals)
    dts[[file]] <- dt
  }
  
  dt <- dts[[1]][, c('cutoff'), with = F]
  dt$FDR <- 0
  for(file in dts){
    dt$FDR <- dt$FDR + file$FDR
  }
  dt$FDR <- dt$FDR / length(dts)
  dt$FDR <- round(dt$FDR * 100)
  write.table(dt, path.out, sep = '\t', row.names = F, col.names = T, quote = T)
  print("Made table!")
  dt
}

Plot <- function(path.to.ph.random, path.to.ph.real, path.out, Xs, cutoffs, dt, x.label = 'energy cutoff, kkal/mol', column, to.file = T){
  Ys <- dt[dt$cutoff %in% cutoffs, ]$FDR
  if(column == 'loop.length'){
    scales.x <- cutoffs
    scales.y <- Ys
  } else{
    scales.x <- c(cutoffs, min(Xs), max(Xs))
    scales.y <- c(min(dt$FDR), Ys, max(dt$FDR))
  }
  secondary.lines.color <- brewer.pal(12, "Set3")[11]
  dots.color <- brewer.pal(12, "Set3")[4]
  p <- ggplot(dt, aes(cutoff, FDR))
    if(column == 'energy'){
      p <- p + geom_line(size = 1.3)
    } else{
      p <- p + geom_smooth(method="lm", formula = y ~ poly(x, 10), se=F, size = 1.3, color = 'black')
    }
    p <- p +
    scale_x_continuous(breaks = scales.x) +
    scale_y_continuous(breaks = scales.y) +
    xlab(x.label) +
    ylab('FDR, %')
  p <- p + theme_linedraw()
  for(x in cutoffs){
    p <- p + 
      geom_point(data = dt[dt$cutoff == x, ], 
                 aes(x = cutoff, y = FDR), col = dots.color, size = 5, pch = 18) 
  }
  p <- p + theme(panel.grid.minor = element_blank(),
                 panel.grid.major.x = element_line(size = 0.7, linetype = 2, 
                                                   color = rep(secondary.lines.color, length(cutoffs) + 1)),
                 panel.grid.major.y = element_line(size = 0.7, linetype = 2, 
                                                   color = rep(secondary.lines.color, length(cutoffs) + 1)))
  p <- p + theme(text = element_text(size = 20),
                 axis.text.x = element_text(size = 15)) +
           theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
  if(column == 'energy'){
    grob <- grobTree(textGrob("spread cutoff = 10.000 nts", x=0.25,  y=0.8, hjust=0,
                              gp=gpar(col="black", fontsize=15)))
    
  } else{
    grob <- grobTree(textGrob("energy cutoff = -15 kcal/mol", x=0.35,  y=0.5, hjust=0,
                              gp=gpar(col="black", fontsize=15)))
  }
  p <- p + annotation_custom(grob)
  if(to.file){
    pdf(path.out, width = 7, height = 5)
    print(p)
    dev.off()
  }
  else {
    print(p)
  }
  print("plotted!")
}

# energy
print("working with energy...")
cutoffs <- c(-15, -20, -25, -30)
Xs <- c(seq(-15, -30, -1), seq(-30, -100, -10))
#dt <- MakeTable(path.to.ph.random, path.to.ph.real, paste0(folder.out, '/energy_cutoff.txt'), Xs, column = "energy")
dt <- fread(paste0(folder.out, '/energy_cutoff.txt'))
Plot(path.to.ph.random, path.to.ph.real, paste0(folder.out, '/energy_cutoff.pdf'), Xs, 
     cutoffs, dt, x.label = 'energy cutoff, kkal/mol', column = "energy", to.file = T)

# Loop lentgh
print("working with loop")
cutoffs <- c(100, 1000, 10000)
Xs <- c(seq(0, 500, 50), seq(500, 10000, 100))
#dt <- MakeTable(path.to.ph.random, path.to.ph.real, paste0(folder.out, '/loop_length_cutoff.txt'), Xs, column = "loop.length")
dt <- fread(paste0(folder.out, '/loop_length_cutoff.txt'))
Plot(path.to.ph.random, path.to.ph.real, paste0(folder.out, '/loop_length_cutoff.pdf'), Xs, 
     cutoffs, dt, x.label = 'spread, nts', column = 'loop.length', to.file = T)


#+ geom_segment(x = min(x), xend = min(x), y = 0, yend = 1)