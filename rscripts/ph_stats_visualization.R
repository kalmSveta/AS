library(data.table)

path.to.ph <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.tsv'
title <- 'filtered'

dt <- fread(path.to.ph)
print('Number of ph = ')
print(length(unique(dt$id)))

png(paste0('pictures/ph_stats/', title, '_handle_length_distr.png'), width = 10, height = 7, res = 300, units = 'in')
hist(dt$al1_length, xlim = c(10, 30), breaks = 100, xlab = 'Handle length, nts', main = '')
dev.off()

dt$length <- dt$panhandle_right_hand - dt$panhandle_left_hand
png(paste0('pictures/ph_stats/', title, '_loop_length_distr.png'), width = 10, height = 7, res = 300, units = 'in')
hist(dt$length, breaks = 50, xlab = 'Loop length, nts', main = '')
dev.off()

png(paste0('pictures/ph_stats/', title, '_energy_distr.png'), width = 10, height = 7, res = 300, units = 'in')
hist(dt$energy, xlim =  c(-40, -15), breaks = 300, xlab = 'delta G, kcal/mol', main = '')
dev.off()
