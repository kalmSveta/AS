library(data.table)
library(ggplot2)
library(stringr)
options(scipen = 999)
path.to.random <- '../python_scripts/folding_pretty_copy/out/FN_rate2/'
path.to.ph.real <- '../python_scripts/folding_pretty_copy/out/hg19_v2/folding/panhandles.tsv'
path.to.comparisons.real <- '../python_scripts/folding_pretty_copy/out/FN_rate/real_data.tsv'
n.files <- 10

##########################
ph.real <- dim(fread(path.to.ph.real))[1]
comparisons.real <- sum(fread(path.to.comparisons.real)$small)


phs.random <- c()
comparisons.random <- c()
for(i in c(1:n.files)){
  ph.random <- dim(fread(paste0(path.to.random, '/random_panhandles', i, '.tsv')))[1]
  comparison.random <- sum(scan(paste0(path.to.random, '/counts_close_', i, '.txt')))
  phs.random <- c(phs.random, ph.random)
  comparisons.random <- c(comparisons.random, comparison.random)
}

normalized.counts.random <- phs.random / comparisons.random
normalized.counts.real <- ph.real / comparisons.real 
fraction <- normalized.counts.random / normalized.counts.real
critical.value <- abs(qt(0.05, df = length(fraction) - 1))
ci <- critical.value * sd(fraction) / sqrt(length(fraction))

print('Condifence interval for FN rate = ')
print(paste0('(', round(100 * (mean(fraction) - ci), 1), ' , ', round(100 * (mean(fraction) + ci), 1), ')%'))

#####################
ph.real <- fread(path.to.ph.real)
ph.real$type <- 'real'
ggplot(ph.real, aes(y = energy)) + 
  geom_boxplot() +
  ylim(c(-100, -15))

files <- list.files(path.to.random, pattern = 'random*', full.names = T)
ph.random <- fread(files[1])
ph.random$type <- 'random'
ggplot(ph.random, aes(y = energy)) + 
  geom_boxplot() +
  ylim(c(-100, -15))

dt <- rbind(ph.real[, c('energy', 'type'), with = F], ph.random[, c('energy', 'type'), with = F])
ggplot(dt, aes(y = energy, x = type)) +
  geom_boxplot() +
  ylim(c(-100, -15))

############### heat map
ph.real <- dim(ph.real[ph.real$energy < -17, ])[1]
phs.random <- c()
for(i in c(1:n.files)){
  ph.random <- fread(paste0(path.to.random, '/random_panhandles', i, '.tsv'))
  ph.random <- dim(ph.random[ph.random$energy < -17, ])[1]
  phs.random <- c(phs.random, ph.random)
}


bins <- fread('../python_scripts/folding_pretty_copy/data/ph_with_bins.tsv')
bins <- bins[, c('interval_chr_start_end_strand', 'bins'), with = F]
ph.real <- fread(path.to.ph.real)
ph.real <- merge(ph.real, bins, by.y = 'interval_chr_start_end_strand', by.x = 'interval1')
ph.real <- merge(ph.real, bins, by.y = 'interval_chr_start_end_strand', by.x = 'interval2')
ph.real[ , `:=`( COUNT = .N) , by = list(ph.real$bins.x, ph.real$bins.y)]
ph.real$length1 <- as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 2])
ph.real$length2 <- as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 2])
ggplot(ph.real, aes(x = bins.x, y = bins.y, fill = COUNT)) +
  geom_tile()
dim(ph.real[ph.real$bins.x != 49 & ph.real$bins.y != 49, ])
dim(ph.real[ph.real$length1 >= 50 & ph.real$length1 <= 1000 & ph.real$length2 >= 50 & ph.real$length2 <= 1000, ])

ph.random <- fread(paste0(path.to.random, '/random_panhandles', i, '.tsv'))
ph.random <- merge(ph.random, bins, by.y = 'interval_chr_start_end_strand', by.x = 'interval1')
ph.random <- merge(ph.random, bins, by.y = 'interval_chr_start_end_strand', by.x = 'interval2')
ph.random[ , `:=`( COUNT = .N) , by = list(ph.random$bins.x, ph.random$bins.y)]
ph.random$length1 <- as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 2])
ph.random$length2 <- as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 2])
ggplot(ph.random, aes(x = bins.x, y = bins.y, fill = COUNT)) +
  geom_tile()
dim(ph.random[ph.random$bins.x != 49 & ph.random$bins.y != 49, ])
dim(ph.random[ph.random$length1 >= 50 & ph.random$length1 <= 1000 & ph.random$length2 >= 50 & ph.random$length2 <= 1000, ])


################################  FDR from length cutoff
ph.real <- fread(path.to.ph.real)
ph.random <- fread('../python_scripts/folding_pretty_copy/out/FN_rate/random_panhandles1.tsv')
print(dim(ph.random)[1] / dim(ph.real)[1])

ph.random$length1 <- as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 2])
ph.random$length2 <- as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 2])

ph.real$length1 <- as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 2])
ph.real$length2 <- as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 3]) - as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 2])

reals <- c()
randoms <- c()
length.cutoffs <- c(10, 11, 12, 13, 15, 17, 20, 30, 50, 100, 300, 500, 1000, 2000, 3000)
for(l in length.cutoffs){
  randoms <- c(randoms, dim(ph.random[ph.random$length1 <= l & ph.random$length2 <= l, ])[1])
  reals <- c(reals, dim(ph.real[ph.real$length1 <= l & ph.real$length2 <= l, ])[1])
}
png('../python_scripts/folding_pretty_copy/out/FDR_from_length_cutoff_upper.png', width = 7, height = 5, units = 'in', res = 300)
plot(x = length.cutoffs, y = randoms/reals)
dev.off()


reals <- c()
randoms <- c()
length.cutoffs <- c(seq(10, 30, 1), seq(31, 100, 5), seq(101, 500, 10), seq(501, 1000, 100), seq(1001, 3000, 1000))
for(l in length.cutoffs){
  randoms <- c(randoms, dim(ph.random[ph.random$length1 >= l & ph.random$length2 <= l, ])[1])
  reals <- c(reals, dim(ph.real[ph.real$length1 >= l & ph.real$length2 <= l, ])[1])
}
png('../python_scripts/folding_pretty_copy/out/FDR_from_length_cutoff_down.png', width = 7, height = 5, units = 'in', res = 300)
plot(x = length.cutoffs, y = randoms/reals)
dev.off()

################################## FDR from energy cutoff
reals <- c()
randoms <- c()
energy.cutoffs <- c(-15, -16, -20, -25, -30, -50, -100, -200, -300, -400)
for(l in energy.cutoffs){
  randoms <- c(randoms, dim(ph.random[ph.random$energy >= l, ])[1])
  reals <- c(reals, dim(ph.real[ph.real$energy >= l, ])[1])
}
png('../python_scripts/folding_pretty_copy/out/FDR_from_energy_cutoff_upper.png', width = 7, height = 5, units = 'in', res = 300)
plot(x = energy.cutoffs, y = randoms/reals)
dev.off()

reals <- c()
randoms <- c()
energy.cutoffs <- c(seq(-15, -130, -5))
#energy.cutoffs <- c(seq(-15, -, -5), seq(-31, -100, -10), seq(-101, -400, -30))
for(l in energy.cutoffs){
  randoms <- c(randoms, dim(ph.random[ph.random$energy <= l, ])[1])
  reals <- c(reals, dim(ph.real[ph.real$energy <= l, ])[1])
}
dt <- data.table('energy cutoff, kkal/mol' = energy.cutoffs, 'FDR' = randoms / reals)
dt$FDR <- round(dt$FDR * 100)
y <- dt[dt$`energy cutoff, kkal/mol` == -30, ]$FDR
p <- ggplot(dt, aes(`energy cutoff, kkal/mol`, FDR)) +
  geom_line() +
  geom_point(data = dt[dt$`energy cutoff, kkal/mol` == -30, ], aes(x = `energy cutoff, kkal/mol`, y = FDR), col = "red", size = 3) +
  scale_x_continuous(breaks = c(-30, -50, -100)) +
  scale_y_continuous(breaks = c(20, y, 40, 60, 80)) + 
  theme(axis.text.x = element_text(color = c("red", "black", "black")),
        axis.ticks.x = element_line(color = c("red", "black", "black"))) +
  theme(axis.text.y = element_text(color = c("black", "red", "black", "black", "black")),
        axis.ticks.y = element_line(color = c("black", "red", "black", "black", "black"))) +
  xlab('energy cutoff, kkal/mol') +
  ylab('FDR, %')

p <- p + theme(panel.grid.minor = element_blank(),
              panel.grid.major.x = element_line(color = c(NA, "white", "white")),
              panel.grid.major.y = element_line(color = c("white", NA, "white", "white", "white")))

p <- p + theme(axis.text.x = element_text(color = "grey20", size = 12),
           axis.text.y = element_text(color = "grey20", size = 12),  
           axis.title.x = element_text(color = "grey20", size = 12),
           axis.title.y = element_text(color = "grey20", size = 12))

png('../python_scripts/folding_pretty_copy/out/FDR_from_energy_cutoff_down.png', width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()

################################## n struct for every pair of lengths
ph.real[ , `:=`( COUNT = .N) , by = list(ph.real$length1, ph.real$length2)]
x <- ph.real[, c('length1', 'length2', 'COUNT'), with = F]
x <- x[!duplicated(x), ]
png('../python_scripts/folding_pretty_copy/out/n_struct_length_real.png', width=7, height=5, res = 300, units = 'in')
plot((x$length1 + x$length2)/2, x$COUNT)
dev.off()

x$length <- (x$length1 + x$length2)/2
png('../python_scripts/folding_pretty_copy/out/hist_n_struct_length_real.png', width=7, height=5, res = 300, units = 'in')
hist(x$length, breaks = 100)
dev.off()

m <- max(x$COUNT)

ph.random[ , `:=`( COUNT = .N) , by = list(ph.random$length1, ph.random$length2)]
x <- ph.random[, c('length1', 'length2', 'COUNT'), with = F]
x <- x[!duplicated(x), ]
png('../python_scripts/folding_pretty_copy/out/n_struct_length_random.png', width=7, height=5, res = 300, units = 'in')
plot((x$length1 + x$length2)/2, x$COUNT, ylim = c(0, m))
dev.off()

x$length <- (x$length1 + x$length2)/2
png('../python_scripts/folding_pretty_copy/out/hist_n_struct_length_random.png', width=7, height=5, res = 300, units = 'in')
hist(x$length, breaks = 100)
dev.off()
############################## n subopt struct for every pair of intrevals
ph.real$inervals <- paste(ph.real$interval1, ph.real$interval2, sep = '_')
x <- as.data.table(t(table(ph.real$inervals)))
png('../python_scripts/folding_pretty_copy/out/suboptimal_length.png', width=7, height=5, res = 300, units = 'in')
x$length1 <- as.numeric(str_split_fixed(x$V2, '_', 4)[, 3]) - as.numeric(str_split_fixed(x$V2, '_', 4)[, 2])
x$length2 <- as.numeric(str_split_fixed(x$V2, '_', 8)[, 7]) - as.numeric(str_split_fixed(x$V2, '_', 8)[, 6])
plot((x$length1 + x$length2)/2, y = log(x$N), ylab = 'log(# suboptimal structures)', xlab = '(interval1+interval2) / 2')
dev.off()

x$length <- (x$length1 + x$length2) / 2
hist(x$length, breaks = 100)

#################################### FDR from loop length
ph.real$interval1_start <- as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 2])
ph.real$interval1_end <- as.numeric(str_split_fixed(ph.real$interval1, '_', 4)[, 3])
ph.real$interval2_start <- as.numeric(str_split_fixed(ph.real$interval2, '_', 4)[, 2])
ph.real$loop_start <- ph.real$interval1_start + ph.real$end_al1
ph.real$loop_end <- ph.real$interval2_start + ph.real$start_al2
ph.real$loop_length <- ph.real$loop_end - ph.real$loop_start

ph.random$interval1_start <- as.numeric(str_split_fixed(ph.random$interval1, '_', 4)[, 2])
ph.random$interval2_start <- as.numeric(str_split_fixed(ph.random$interval2, '_', 4)[, 2])
ph.random$loop_start <- ph.random$interval1_start + ph.random$end_al1
ph.random$loop_end <- ph.random$interval2_start + ph.random$start_al2
ph.random$loop_length <- ph.random$loop_end - ph.random$loop_start

library(hexbin)
library(RColorBrewer)

reals <- c()
randoms <- c()
energy.cutoffs <- c(seq(0, 500, 50), seq(500, 10000, 100))
for(l in energy.cutoffs){
  randoms <- c(randoms, dim(ph.random[ph.random$loop_length <= l, ])[1])
  reals <- c(reals, dim(ph.real[ph.real$loop_length <= l, ])[1])
}
dt <- data.table('loop cutoff, nts' = energy.cutoffs, 'FDR' = randoms / reals)
dt$FDR <- round(dt$FDR * 100)

p <- ggplot(dt, aes(`loop cutoff, nts`, FDR)) +
  geom_line() +
  xlab('loop cutoff, nts') +
  ylab('FDR, %')



png('../python_scripts/folding_pretty_copy/out/FDR_from_energy_cutoff_down.png', width = 7, height = 5, units = 'in', res = 300)
print(p)
dev.off()

x <- ph.real[ph.real$energy > -60, ]$loop_length
y <- ph.real[ph.real$energy > -60, ]$energy
bin<-hexbin(x, y, xbins=100)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))

png("../python_scripts/folding_pretty_copy/out/FDR_figures/energy_from_loop_length.png", width = 7, height = 7, units = 'in', res = 300)
plot(bin, main="" , colramp=my_colors , legend=F, xlab = 'loop length, nts', ylab = 'energy, kkal/mol')
dev.off()

png("../python_scripts/folding_pretty_copy/out/FDR_figures/hist_distance_between_conins.png", width = 7, height = 5, units = 'in', res = 300)
hist(ph.real$interval2_start - ph.real$interval1_end, breaks = 100, main = '', xlab = 'Distance between conserved intronic intervals, nts')
dev.off()

png("../python_scripts/folding_pretty_copy/out/FDR_figures/hist_loop_length.png", width = 7, height = 5, units = 'in', res = 300)
hist(ph.real$loop_length, breaks = 100, main = '', xlab='loop length')
dev.off()

png("../python_scripts/folding_pretty_copy/out/FDR_figures/hist_loop_length_energy_cutoff_30.png", width = 7, height = 5, units = 'in', res = 300)
hist(ph.real[ph.real$energy <= -30, ]$loop_length, breaks = 100, main = '', xlab='loop length')
dev.off()


###################
shuffled.conins <- fread('../python_scripts/folding_pretty_copy/out/intervals_shuffled_1.tsv')
shuffled.conins <- shuffled.conins[, c('interval_chr_start_end_strand', 'sequences'), with = F]
shuffled.ph <- fread('../python_scripts/folding_pretty_copy/out/random_panhandles1.tsv')
merged <- merge(shuffled.ph, shuffled.conins, by.x = 'interval1', by.y = 'interval_chr_start_end_strand')
colnames(merged)[colnames(merged) == 'sequences'] <- 'interval1_sequence'
merged <- merge(merged, shuffled.conins, by.x = 'interval2', by.y = 'interval_chr_start_end_strand')
colnames(merged)[colnames(merged) == 'sequences'] <- 'interval2_sequence'

