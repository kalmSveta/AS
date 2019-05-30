library(data.table)
options(scipen = 999)
path.to.ph <- '../python_scripts/folding_pretty_copy/out/folding/panhandles_preprocessed.tsv'
mut <- fread('../python_scripts/compensatory_copy/out/not_filtered/ph_mutation_amount.tsv')


ph <- fread(path.to.ph)
MutationalScore <- function(ph){
  mut$mutational_score <- mut$n_mut * (-1) + mut$n_comp * (+10)
  ph <- merge(ph, mut[, c('ph_id', 'mutational_score'), with = F], by.x = 'id', by.y = 'ph_id', all.x = T)
  ph[is.na(ph$mutational_score), ]$mutational_score <- 0
  ph
}
ph <- MutationalScore(ph)
write.table(ph, path.to.ph, sep = '\t', row.names = F, quote = F)
