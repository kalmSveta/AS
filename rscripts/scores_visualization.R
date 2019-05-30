library(data.table)
library(ggplot2)
dt <- fread('../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed.tsv')
filtered <- fread('../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed_filtered.tsv')
dt <- dt[dt$id %in% filtered$id, ]
known.ids <- c(127812, 70913, 579130, 106707,107608, 107609)
names(known.ids) <- c('SF1', 'ENAH', 'DNM1', 'ATE1_1', 'ATE1_2', 'ATE1_3')

dt[dt$gene_name == 'GFRA1',]$gene_name <- 'ATE1_2'

for(score in c("SplicingScore", "mutational_score", "conservation_score", "kmerScore", "Energy_and_freq_Score")){
  svg(paste0('~/Desktop/Thesis_pictures/Results/', score, '.svg'), width = 7, height = 7)
  known.scores <- dt[dt$id %in% known.ids, ][[score]]
  names(known.scores) <- dt[dt$id %in% known.ids, ]$gene_name
  hist(dt[[score]], xlab = 'score', main = score)
  abline(v = known.scores, col = 'red')
  text(known.scores, 60000, names(known.scores), cex = 0.6, pos = 1, col = "red") 
  dev.off()
}
