library(stringr)
revcom <- function(s){
  s <- chartr("ATGC","TACG",s)
  splits <- strsplit(s, "")[[1]]
  reversed <- rev(splits)
  final_result <- paste(reversed, collapse = "")
  final_result
}

revSTRUCT <- function(s){
  left <- substr(s, 1, str_locate(s, '\\)')[1, 'start'] - 1)
  right <- substr(s, str_locate(s, '\\)')[1, 'start'], nchar(s))
  paste0(gsub('\\)', '\\(', right), gsub('\\(', '\\)', left))
}

dt <- fread('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv')
dt$revc1 <- unlist(lapply(dt$alignment1, function(x) revcom(x)))
dt$revc2 <- unlist(lapply(dt$alignment2, function(x) revcom(x)))
dt[dt$strand == '-', ]$alignment1 <- dt[dt$strand == '-', ]$revc1
dt[dt$strand == '-', ]$alignment2 <- dt[dt$strand == '-', ]$revc2



dt$rev <- unlist(lapply(dt$structure, function(x) revSTRUCT(x)))
dt[dt$strand == '-', ]$structure <- dt[dt$strand == '-', ]$rev
dt <- dt[, !names(dt) %in% c('rev'), with = F]
write.table(dt, '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered_struct.tsv', 
            sep = '\t', row.names = F, col.names = T, quote = F)
