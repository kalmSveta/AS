library(data.table)
library(R.utils)
library(stringr)

path.to.mut.files <- ('../python_scripts/compensatory_copy/out/handles_and_mut/filtered/')
path.to.ph <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.tsv'



files <- list.files(path.to.mut.files, '*txt.gz', full.names = T)

mut.list <- list()
for(file in files){
  print(file)
  mut.dt <- fread(file)
  mut.dt <- mut.dt[, c(1, 2, 3, 4, 5, 7, 8, 9, 10), ]
  mut.list[[file]] <- mut.dt
}

mut.dt <- do.call(rbind, mut.list)
mut.dt$V1 <- paste('chr', mut.dt$V1, sep = '')
mut.dt$panhandle_id <- as.numeric(str_split_fixed(mut.dt$V5, '_', 2)[, 1])
mut.dt$panhandle_hand <- str_split_fixed(mut.dt$V5, '_', 2)[, 2]
mut.dt <- mut.dt[, c('V1', 'panhandle_id', 'panhandle_hand', 'V7', 'V8', 'V9', 'V10'), with = F]
colnames(mut.dt) <- c('chr', 'panhandle_id', 'panhandle_hand', 'mut_coord', 'mutation_info', 'mut_from', 'mut_to')

ph.dt <- fread(path.to.ph)
mut.dt <- merge(mut.dt, ph.dt, by.x = c('chr', 'panhandle_id'), by.y = c('chr', 'id'), all.x = T)
mut.dt <- mut.dt[!is.na(mut.dt$gene_id), ]
mut.dt <- mut.dt[grep(',', mut.dt$mut_to, invert = T), ]
write.table(mut.dt, './filtered_panhandles_with_populational_mutations.tsv', sep = '\t', row.names=F, quote=F)
