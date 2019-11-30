library(data.table)
library(stringr)
options(scipen = 999)
anno <- '../conservative_features/gencode.v19.annotation.gtf'


anno.dt <- fread(anno, skip = 5, header = F)
anno.dt$gene.id <- str_split_fixed(anno.dt$V9, '; ', 2)[, 1]
# select coding genes
coding.genes <- unique(anno.dt[anno.dt$V3 == 'CDS', ]$gene.id)
anno.dt <- anno.dt[anno.dt$V3 == 'gene' & anno.dt$gene.id %in% coding.genes, ]
# select not intersected genes
anno.dt <- anno.dt[, c('V1', 'V4', 'V5','gene.id', 'V6', 'V7'), with = F]
anno.dt$gene.id <- gsub('gene_id ', '', anno.dt$gene.id)
anno.dt$gene.id <- gsub('\\"', '', anno.dt$gene.id)
write.table(anno.dt, 'tmp.txt', sep = '\t', row.names = F, col.names = F, quote = F)
dt <- read.delim(pipe("bedtools intersect -a tmp.txt -b tmp.txt -s -c"), header = F)
dt <- dt[dt$V7 == 1, ]
dt <- dt[, !names(dt) %in% c('V7')]
write.table(dt, '../conservative_features/not_intersected_coding_genes.bed', sep = '\t', row.names = F, col.names = F, quote = F)
