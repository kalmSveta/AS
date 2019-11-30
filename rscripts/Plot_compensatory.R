#!/usr/bin/Rscript
library(dplyr)
library(ggpubr)
options(scipen = 999)

comp.shuffled <- '../python_scripts/compensatory_copy/out/N_compensatory_shuffled.txt'
pairs.shuffled <- '../python_scripts/compensatory_copy/out/N_pairs_with_SNP_shuffled.txt'
comp <- 40
pairs <- 82994
title <- '100 iterations'

# Plot_compensatory.R N_pairs_with_SNP_shuffled N_compensatory_shuffled N_pairs_with_SNP N_compensatory title
args <- commandArgs(trailingOnly = TRUE)
comp.shuffled <- args[1]
pairs.shuffled <- args[2]
comp <- as.numeric(args[3])
pairs <- as.numeric(args[4])
title <- args[5]

comp.shuffled <- as.numeric(scan(comp.shuffled, what = 'numeric', quote = ""))
pairs.shuffled <- as.numeric(scan(pairs.shuffled, what = 'numeric'))
real.ratio <- comp / pairs * 100

Zscore <- (real.ratio - mean(comp.shuffled / pairs.shuffled * 100)) / sd(comp.shuffled / pairs.shuffled * 100)
pvalue <- pnorm(q = Zscore, lower.tail = F)
pdf(paste0('Distribution_of_compensatory_shuffled_in_', title, '.pdf'))
width.min <- min((comp.shuffled / pairs.shuffled * 100))
width.max <- max((comp.shuffled / pairs.shuffled * 100))
h <- hist(comp.shuffled / pairs.shuffled * 100, breaks = 30, xlim = c(width.min - 0.005, max(width.max, real.ratio) + 0.005), 
     main = paste0('Percantage of panhandles bp with compensatory mutations \namong panhandles bp with mutations \n', 
                   title, ', Zscore = ', round(Zscore, 2), ', p-value = ', round(pvalue, 2)), 
     xlab = '% bp with compensatory mutations')
height <- max(h$counts)
abline(v = real.ratio, col = 'red')
text(real.ratio * 0.7, height * 0.8,  cex = 1.3, pos = 4, "Observed", col = 'red')
text(width.min + 0.003, height * 0.8, cex = 1.3, pos = 4, "Expected", col = 'black')
dev.off()

# 
# shapiro.test(comp.shuffled / pairs.shuffled * 100)
# png('qqplot.png', width = 5, height = 5, units = 'in', res = 300)
# ggqqplot(comp.shuffled / pairs.shuffled)
# dev.off()


