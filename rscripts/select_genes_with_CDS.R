#!/usr/bin/Rscript
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

# Usage: Rscript select_exons_with_CDS.R path.to.anno out.path

path.to.anno <- '../conservative_features/gencode.v19.annotation.gtf'
out.path <- '../conservative_features/coding_genes.bed'

path.to.anno <-  args[1]
out.path <- args[2]


SelectCodingExons <- function(path.to.anno, out.path){
  dt <- fread(path.to.anno, skip = 5)
  dt$gene_name <- str_split_fixed(dt$V9, '; ', 6)[, 5]
  dt$gene_id <- str_split_fixed(dt$V9, '; ', 2)[, 1]  
  CDS <- dt[dt$V3 == 'CDS', ] 
  dt <- dt[dt$V3 == 'gene',]   
  dt <- dt[dt$gene_id %in% CDS$gene_id, ]
  dt <- dt[, c("V1", "V4", "V5", 'gene_id', 'gene_name', 'V7'), with = F]
  write.table(dt, out.path, sep = '\t', col.names = F, row.names = F, quote = F)  
}

SelectCodingExons(path.to.anno, out.path)