.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(stringr)
options(scipen = 999)


energy.cutoffs <- c(-15, -20, -25, -30)
energy.colors <- brewer.pal(n = 12, "Paired")[c(4, 7, 8, 6)]
names(energy.colors) <- paste('<=', energy.cutoffs, sep = '')


# path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/panhandles_preprocessed_filtered.tsv'
# path.out <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/stats/'
# title <- 'filtered'
# path.to.conin <- '../python_scripts/folding_pretty_copy/out/hg19_strand_specific/intervals_with_seqs.tsv'

path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv'
path.out <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/stats/'
dir.create(file.path(path.out), showWarnings = FALSE)
title <- 'filtered'
path.to.conin <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/intervals_with_seqs.tsv'

conin.dt <- fread(path.to.conin)
dt <- fread(path.to.ph)
print('Number of ph = ')
print(length(unique(dt$id)))

# handle length (suppl)
pdf(paste0(path.out, title, '_handle_length_distr.pdf'), width = 10, height = 7)
ggplot(data = dt, aes(al1_length)) + 
  geom_histogram(bins = 20, fill = adjustcolor(brewer.pal(12, "Set3")[6], alpha.f = 0.4), color = brewer.pal(12, "Set3")[4]) + 
  xlim(10, 30) + 
  theme_linedraw() +
  xlab('Handle length, nts') +
  ylab('count') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) 
dev.off()

# spread
dt$loop.length <- dt$panhandle_right_hand - dt$panhandle_left_hand
pdf(paste0(path.out, title, '_loop_length_distr.pdf'), width = 10, height = 7)
ggplot(data = dt, aes(loop.length)) + 
  geom_histogram(bins = 30, fill = adjustcolor(brewer.pal(12, "Set3")[6], alpha.f = 0.4), color = brewer.pal(12, "Set3")[4]) + 
  xlim(0, 10000) + 
  theme_linedraw() +
  xlab('Spread, nts') +
  ylab('count') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  geom_abline(slope = 0, intercept = 16000, linetype = 'dashed')
dev.off()

#dt$log.energy <- log10(-dt$energy)
cairo_pdf(paste0(path.out, title, '_energy_distr.pdf'), width = 10, height = 7)
ggplot(data = dt, aes(energy)) + 
  geom_histogram(bins = 20, fill = adjustcolor(brewer.pal(12, "Set3")[6], alpha.f = 0.4), color = brewer.pal(12, "Set3")[4]) + 
  theme_linedraw() +
  #xlab(expression(paste('log10(-', Delta, 'G), kcal/mol'))) +
  xlab('- \u0394G, kcal/mol') +
  ylab('frequency') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  scale_x_reverse(limits = c(max(dt$energy), -30))
dev.off()

dt$position <- (dt$panhandle_start - dt$start_gene) / 
  ((dt$end_gene - dt$start_gene) - (dt$panhandle_end - dt$panhandle_start) + 1)
dt[dt$strand == '-', ]$position <- 1- dt[dt$strand == '-', ]$position 
pdf(paste0(path.out, title, '_ph_position_in_gene.pdf'), width = 10, height = 7)
ggplot(data = dt, aes(position)) + 
  geom_histogram(bins = 30, fill = adjustcolor(brewer.pal(12, "Set3")[6], alpha.f = 0.4), 
                 color = brewer.pal(12, "Set3")[4]) + 
  theme_linedraw() +
  xlab('relative position in gene') +
  ylab('count') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) 
dev.off()


dt1 <- dt[, c('start_gene', 'end_gene', 'panhandle_start', 'panhandle_left_hand', 'strand'), with = F]
dt1$position <- (dt1$panhandle_start - dt1$start_gene) / 
  ((dt1$end_gene - dt1$start_gene) - (dt1$panhandle_left_hand - dt1$panhandle_start) + 1)
dt1[dt1$strand == '-', ]$position <- 1- dt1[dt1$strand == '-', ]$position 
dt1 <- dt1[, c('start_gene', 'end_gene', 'position'), with = F]
dt1$handle <- 'left'
dt2 <- dt[, c('start_gene', 'end_gene', 'panhandle_right_hand', 'panhandle_end', 'strand'), with = F]
dt2$position <- (dt2$panhandle_right_hand - dt2$start_gene) / 
  ((dt2$end_gene - dt2$start_gene) - (dt2$panhandle_end - dt2$panhandle_right_hand) + 1)
dt2[dt2$strand == '-', ]$position <- 1- dt2[dt2$strand == '-', ]$position 
dt2 <- dt2[, c('start_gene', 'end_gene', 'position'), with = F]
dt2$handle <- 'right'
dt <- rbind(dt1, dt2)

pdf(paste0(path.out, title, '_handle_position_in_gene.pdf'), width = 10, height = 7)
ggplot(data = dt, aes(position, fill = handle, color = handle)) + 
  geom_histogram(bins = 30, position = 'identity') + 
  theme_linedraw() +
  xlab('relative position in gene') +
  ylab('count') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = c(adjustcolor(brewer.pal(12, "Set3")[6], alpha.f = 0.3), 
                               adjustcolor(brewer.pal(12, "Set3")[5], alpha.f = 0.3))) +
  scale_color_manual(values = c(brewer.pal(12, "Set3")[4], 
                               brewer.pal(12, "Set3")[5]))
dev.off()




## spread and intervals
conin.dt <- conin.dt[, c('start_interval', 'end_interval', 'gene_chr_start_end_strand', 'interval_chr_start_end_strand'), with = F]
interval.distances <- unlist(lapply(unique(conin.dt$gene_chr_start_end_strand), function(gene){
  tmp <- conin.dt[conin.dt$gene_chr_start_end_strand == gene, ]
  tmp <- merge(tmp, tmp, by = 'gene_chr_start_end_strand', allow.cartesian=TRUE)
  tmp <- tmp[tmp$interval_chr_start_end_strand.y != tmp$interval_chr_start_end_strand.x, ]
  tmp <- tmp[tmp$start_interval.y >= tmp$end_interval.x, ]
  tmp$start_interval.y - tmp$end_interval.x
}))
interval.distances <- interval.distances[interval.distances <= 10000]
write(interval.distances, 'interal_distances.txt', sep = '\t')
interval.distances <- scan('interal_distances.txt')

# ratio of densities
bins <- 15
dt$loop.length <- dt$panhandle_right_hand - dt$panhandle_left_hand
spread.density <- as.data.table(density(dt$loop.length, from = min(dt$loop.length), to = max(dt$loop.length), n = bins)[c('x', 'y')])
interval.distances.density <- as.data.table(density(interval.distances, from = min(dt$loop.length), to = max(dt$loop.length), n = bins)[c('x', 'y')])
merged.density <- merge(interval.distances.density, spread.density, by = 'x', suffixes = c('.interval', '.spread'))
merged.density$ratio <- merged.density$y.interval / merged.density$y.spread
breaks = c(seq(0, 1.2, by=0.4), 1)
labels = as.character(breaks)
pdf(paste0(path.out, title, '_density_ratio_spread_and_interval_distance.pdf'))
ggplot(merged.density, aes(x, ratio)) +
  geom_line(size = 1, color = brewer.pal(12, "Set3")[7]) +
  xlab('spread, nts') +
  theme_linedraw() +
  geom_abline(slope = 0, intercept = 1, linetype = 'dashed', color = 'red') + 
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) +
  scale_y_continuous(breaks = breaks, labels = labels,
                     name = "interval distance density / spread density") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# two densities one plot
interval.distances <- scan('interal_distances.txt')
dt$spread <- dt$panhandle_right_hand - dt$panhandle_left_hand
frequency.dt <- data.table(distances = interval.distances, type = 'distance between conserved intronic intervals')
spread <- data.table(distances = dt$spread, type = 'spread')
frequency.dt <- rbind(frequency.dt, spread)
frequency.dt$type <- factor(frequency.dt$type)

pdf(paste0(path.out, title, '_spread_hist.pdf'), width = 7, height = 5)
ggplot(frequency.dt, aes(distances, fill = type, color = type)) +
  #geom_density(position = 'identity', adjust = 0.1, alpha = 0.3) +
  geom_histogram(aes(y=..density..), position = 'identity', alpha = 0.3,  bins = 50) +
  #geom_histogram(position = 'identity', alpha = 0.1,  bins = 20) +
  scale_fill_manual(name = '', values = brewer.pal(12, "Set3")[c(5, 6)]) +
  scale_color_manual(name = '', values = brewer.pal(12, "Set3")[c(5, 6)]) +
  theme_linedraw() +
  xlab('nts') +
  ylab('density') +
  xlim(min(frequency.dt$distances), 10000) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = c(0.55,0.85), legend.direction = "vertical")
dev.off()



# panhandle position in gene and conin
min.conins <- aggregate(start_interval ~ gene_chr_start_end_strand, conin.dt, min)
max.conins <- aggregate(end_interval ~ gene_chr_start_end_strand, conin.dt, max)
range.conins <- merge(min.conins, max.conins, by = 'gene_chr_start_end_strand')

dt1 <- dt[, c('start_gene', 'end_gene', 'panhandle_start', 'panhandle_end', 'strand', 'chr', 'energy', 'id'), with = F]
dt1$gene_chr_start_end_strand <- paste(dt1$chr, dt1$start_gene, dt1$end_gene, dt1$strand, sep = '_')
dt1 <- merge(dt1, range.conins, by = 'gene_chr_start_end_strand')

MakeCategories2 <- function(dt, what = 'energy', categories){
  dt.new <- copy(dt)
  dt.new[, paste0(what, '.bin') := categories[1]]
  for(cat in categories[2:length(categories)]){
    dt.tmp <- dt.new[dt.new[[what]] <= cat, ]
    dt.tmp[, paste0(what, '.bin') := cat]
    dt.new <- rbind(dt.new, dt.tmp)
  }
  #dt.new <- dt.new[!duplicated(dt.new), ]
  dt.new
}

dt1$position <- (dt1$panhandle_start - dt1$start_interval) / 
  ((dt1$end_interval - dt1$start_interval) - (dt1$panhandle_end - dt1$panhandle_start) + 1)
dt1[dt1$strand == '-', ]$position <- 1 - dt1[dt1$strand == '-', ]$position 
dt1 <- MakeCategories2(dt1, what = 'energy', energy.cutoffs)
dt1$energy.bin <- paste('<=', dt1$energy.bin, sep = '')
dt1$energy.bin <- factor(dt1$energy.bin)
cairo_pdf(paste0(path.out, title, '_panhandle_position_in_gene_interval.pdf'), width = 10, height = 7)
ggplot(data = dt1, aes(position, fill = energy.bin, color = energy.bin)) + 
  geom_density(position = 'identity', adjust = 3, alpha = 0.05) +
  geom_hline(yintercept=0, colour="white", size=1) + 
  geom_vline(xintercept=0, colour="white", size=1) +
  geom_vline(xintercept=1, colour="white", size=1) +
  theme_linedraw() +
  #coord_cartesian(xlim = c(0.045, 0.953)) +
  scale_fill_manual(name = '\u0394G, kcal/mol', values = energy.colors) +
  scale_color_manual(name = '\u0394G, kcal/mol', values = energy.colors) +
  xlab('relative position in gene (from the first interval to the last interval)') +
  ylab('density') + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15))  
dev.off()
