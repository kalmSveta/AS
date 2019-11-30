#!/usr/bin/Rscript
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

# Usage: Rscript select_exons_with_CDS.R path.to.anno out.path

path.to.anno <- '../conservative_features/gencode.v19.annotation.gtf'
out.path <- '../conservative_features/coding_exons.bed'

path.to.anno <-  args[1]
out.path <- args[2]



SelectCodingExons <- function(path.to.anno, out.path){
  dt <- fread(path.to.anno, skip = 5)
  dt$coding <- F
  dt$id <- c(1:dim(dt)[1])
  CDS.ids <- dt[dt$V3 == 'CDS', ]$id
  CDS.ids <- CDS.ids - 1
  dt[dt$id %in% CDS.ids, ]$coding <- T
  dt <- dt[dt$coding, ]
  dt <- dt[, c('V1', 'V4', 'V5', 'id', 'V6', 'V7'), with = F]
  dt$V6 <- 1
  dt$V6 <- as.numeric(dt$V6)
  write.table(dt, out.path, sep = '\t', col.names = F, row.names = F, quote = F)  
}

SelectCodingExons(path.to.anno, out.path)