#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(data.table)
library(R.utils)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
path.to.mut.folder <- args[1]
path.to.ph <- args[2]

#path.to.mut.folder <- '../1000genomes/new/hg19_ss_flanks/phs/'
#path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv'


list.of.files <- file.info(list.files(path.to.mut.folder, full.names = T))
list.of.non.empty.files <- subset(list.of.files, size != 0)
files <- rownames(list.of.non.empty.files)

mut.list <- list()
for(file in files){
  print(file)
  mut.dt <- fread(file)
  mut.dt <- mut.dt[, c(1:11), with = F]
  mut.list[[file]] <- mut.dt
}

mut.dt <- do.call(rbind, mut.list)
mut.dt$V1 <- paste('chr', mut.dt$V1, sep = '')
mut.dt$panhandle_id <- str_split_fixed(mut.dt$V4, '_', 2)[, 1]
mut.dt$panhandle_hand <- str_split_fixed(mut.dt$V4, '_', 2)[, 2]
mut.dt <- mut.dt[, c('V1', 'panhandle_id', 'panhandle_hand', 'V8', 'V9', 'V10', 'V11'), with = F]
colnames(mut.dt) <- c('chr', 'panhandle_id', 'panhandle_hand', 'mut_coord', 'mutation_info', 'mut_from', 'mut_to')

ph.dt <- fread(path.to.ph)
ph.dt$id <- as.character(ph.dt$id)
mut.dt <- merge(mut.dt, ph.dt[, c(1:17), with = F], by.x = c('chr', 'panhandle_id'), by.y = c('chr', 'id'), all.x = T)
#mut.dt <- mut.dt[!is.na(mut.dt$gene_id), ]
mut.dt <- mut.dt[grep(',', mut.dt$mut_to, invert = T), ]
mut.dt$mutation_info <- c(1:nrow(mut.dt))
write.table(mut.dt, gsub('.tsv', '_with_populational_mutations.tsv', path.to.ph), sep = '\t', row.names=F, quote=F)
