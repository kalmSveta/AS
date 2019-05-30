library(data.table)
library(ggplot2)

path.to.random <- '../python_scripts/folding_pretty_copy/out/FN_rate2/'
path.to.ph.real <- '../python_scripts/folding_pretty_copy/out/folding/panhandles.tsv'
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

#####################33
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

###############
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
