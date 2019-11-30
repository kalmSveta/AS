library(data.table)
library(stringr)
library(mltools)
options(scipen = 999)


SelectWorkingSample <- function(dt, genes){
  new.dt.name <- gsub('.bed', '_sample.bed', dt)
  system(paste0('bedtools intersect -a ', dt, ' -b ', genes, ' -f 0.9 -u -s > ', new.dt.name))
  return(new.dt.name)
}

# phs start-end
path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
ph <- fread(paste0(path.to.ph, '/panhandles_preprocessed_filtered_dGcutoff-15.bed'))
ph <- ph[, c(1:6), with = F]
write.table(ph, paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6'), 
            sep = '\t', row.names = F, col.names = F, quote = F)
system(paste0('bedtools sort -i ', paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6'), ' > tmp'))
system(paste0('mv tmp ', paste0(path.to.ph, '/panhandles_preprocessed_filtered.bed6')))

# phs left.handle-right.handle
path.to.ph <- '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/'
ph <- fread(paste0(path.to.ph, '/panhandles_preprocessed_filtered.tsv'))
ph$score <- 0
ph$name <- paste(ph$id, ph$energy, sep = '_')
ph <- ph[, c('chr', 'panhandle_left_hand', 'panhandle_right_hand', 'name', 'score', 'strand'), with = F]
write.table(ph, paste0(path.to.ph, '/panhandles_preprocessed_filtered_inner.bed6'), 
            sep = '\t', row.names = F, col.names = F, quote = F)
system(paste0('bedtools sort -i ', paste0(path.to.ph, '/panhandles_preprocessed_filtered_inner.bed6'), ' > tmp'))
system(paste0('mv tmp ', paste0(path.to.ph, '/panhandles_preprocessed_filtered_inner.bed6')))

# handles
MakeHandlesBed <- function(dt){
  dt.left <- dt[, c('chr', 'panhandle_start', 'panhandle_left_hand', 'id', 'strand'), with = F]
  dt.right <- dt[, c('chr', 'panhandle_right_hand', 'panhandle_end', 'id', 'strand'), with = F]
  colnames(dt.left) <- c('chr', 'start', 'stop', 'name', 'strand')
  colnames(dt.right) <- c('chr', 'start', 'stop', 'name', 'strand')
  dt.left$name <- paste0(dt.left$name, '_left')
  dt.right$name <- paste0(dt.right$name, '_right')
  dt.both <- rbind(dt.left, dt.right)
  dt.both$score <- 1
  dt.both <- dt.both[, c('chr', 'start', 'stop', 'name', 'score', 'strand'), with = F]
  dt.both
}
ph <- fread(paste0(path.to.ph, '/panhandles_preprocessed_filtered.tsv'))
handles <- MakeHandlesBed(ph)
write.table(handles, paste0(path.to.ph, '/panhandles_preprocessed_filtered_handles.bed6'), 
            sep = '\t', row.names = F, col.names = F, quote = F)


# introns simple
system('bedtools sort -i ../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed > tmp')
system('mv tmp ../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed')

# not intersected introns simple
dt <- read.delim(pipe("bedtools intersect -a ../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed -b ../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed -s -c -wa"), header = F)
dt <- dt[dt$V7 == 1, ]
write.table(dt, '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_not_intersected.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)

# short introns
dt <- fread('../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed')
dt <- subset(dt, V3 - V2 <= 10000)
write.table(dt, '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_shorter10000.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)

# introns with and without ph of same length 
dt <- SelectWorkingSample('../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed', 
                          '../conservative_features/not_intersected_coding_genes.bed')
dt <- fread(dt)
introns.with.ph <- read.delim(pipe(paste0('bedtools intersect -a ', 
                                          '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python.bed', 
                                          ' -b ', '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks//panhandles_preprocessed_filtered.bed6',
                                          ' -s -wa -wb -F 01')), 
                              header = F)

introns.with.ph <- as.data.table(introns.with.ph)
introns.with.ph$id <- paste(introns.with.ph$V1, introns.with.ph$V2, introns.with.ph$V3, introns.with.ph$V6, sep = '_')
dt$id <- paste(dt$V1, dt$V2, dt$V3, dt$V6, sep = '_')
dt$ph <- F
dt[dt$id %in% unique(introns.with.ph$id), ]$ph <- T
dt$intron.length <- dt$V3 - dt$V2
dt <- dt[order(dt$intron.length), ]
dt <- dt[!duplicated(dt$id), ]
dt$bin <- bin_data(x = dt$intron.length, bins = 500)

x <- table(dt$bin, dt$ph)
sum(x[, 'FALSE'] < x[, 'TRUE'])
x[x[, 'FALSE'] < x[, 'TRUE'], ]
intron.list <- lapply(rownames(x), function(bin.){
  take <- min(x[bin., ])
  tmp <- dt[dt$bin == bin., ]
  tmp_ph <- subset(tmp, ph)
  tmp_no_ph <- subset(tmp, !ph)
  take_ph <- sample(nrow(tmp_ph), size = take)
  tmp_ph <- tmp_ph[take_ph, ]
  take_no_ph <- sample(nrow(tmp_no_ph), size = take)
  tmp_no_ph <- tmp_no_ph[take_no_ph, ]
  tmp <- rbind(tmp_ph, tmp_no_ph)
  tmp
})
intron.population <- Reduce(rbind, intron.list)
intron.population <- intron.population[, c(1:6), with = F]
write.table(intron.population, '../python_scripts/folding_pretty_copy/data/hg19/noCDS_flanks/introns_python_balanced_set_for_tables.bed',
            sep = '\t', row.names = F, col.names = F, quote = F)




# introns with counts ICGC
counts.dt <- fread('../bec08a02-dd30-11e4-911e-ff7fc254ba06.A06.ssj.tsv')
counts.dt[is.na(counts.dt$V6), ]$V6 <- ''
counts.dt <- counts.dt[counts.dt$V6 == 'GTAG', ]
counts.dt$chr <- str_split_fixed(counts.dt$V1, '_', 4)[, 1]
counts.dt$start <- as.numeric(str_split_fixed(counts.dt$V1, '_', 4)[, 2])
counts.dt$end <- as.numeric(str_split_fixed(counts.dt$V1, '_', 4)[, 3])
counts.dt$strand <- str_split_fixed(counts.dt$V1, '_', 4)[, 4]
counts.dt <- counts.dt[, c('chr', 'start', 'end', 'strand', 'V2'), with = F]
counts.dt <- counts.dt[order(c(counts.dt$chr, counts.dt$start, counts.dt$end)), ]
counts.dt <- counts.dt[!is.na(counts.dt$chr), ]
introns <- copy(counts.dt)
colnames(introns)[ncol(introns)] <- 'counts'
introns$name <- c(1:nrow(introns))
introns <- introns[, c('chr', 'start', 'end', 'name', 'counts', 'strand'), with = F]
write.table(introns, '../intron_counts.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# introns with counts ENCODE
# counts.dt <- fread('../introns/rnaseq_i.bed')
# counts.dt$name <- c(1:nrow(counts.dt))
# counts.dt <- counts.dt[, c('V1', 'V2', 'V3', 'name', 'V7', 'V6'), with = F]
# write.table(counts.dt, '../introns/rnaseq_i_pretty.bed', sep = '\t', row.names = F, col.names = F, quote = F)
files.K562 <- c('ENCFF636QII.A07.gff', 'ENCFF486POD.A07.gff')
files.HepG2 <- c('ENCFF670LIE.A07.gff', 'ENCFF074BOV.A07.gff')
path.to.folder <- '../introns/'
ParseFiles <- function(path.to.folder, file){
  dt <- fread(paste0(path.to.folder, file))
  dt <- subset(dt, V3 == 'intron')
  dt <- subset(dt, grepl('chr', dt$V1))
  dt$inc <- str_split_fixed(dt$V9, '; ', 7)[, 6]
  dt$inc <- gsub('"', '', dt$inc)
  dt$inc <- gsub('nDA', '', dt$inc)
  dt$inc <- gsub(' ', '', dt$inc)
  dt$inc <- as.numeric(dt$inc)
  dt <- dt[, c('V1', 'V4', 'V5', 'inc', 'V7'), with = F]
  colnames(dt) <- c('chr', 'start', 'stop', 'score', 'strand')
  dt
}

dt.list <- lapply(files.K562, function(file) ParseFiles(path.to.folder, file))
dt <- Reduce(function(...) merge(..., by = c('chr', 'start', 'stop', 'strand'), all = T), dt.list)
dt[is.na(dt$score.x), ]$score.x <- 0 
dt[is.na(dt$score.y), ]$score.y <- 0 
dt$score <- dt$score.x + dt$score.y
dt$name <- c(1:nrow(dt))
dt <- dt[, c('chr', 'start', 'stop', 'name', 'score', 'strand'), with = F]
dt <- subset(dt, score >= 10)
write.table(dt, paste0(path.to.folder, 'introns_K562.bed'), sep = '\t', row.names = F, col.names = F, quote = F)

# introns ENCODE short
dt <- fread(paste0(path.to.folder, 'introns_K562.bed'))
dt <- subset(dt, V3 - V2 <= 10000)
write.table(dt, paste0(path.to.folder, 'introns_K562_shorter10000.bed'), 
            sep = '\t', row.names = F, col.names = F, quote = F)


# exons with counts
# counts.dt <- fread('../bec08a02-dd30-11e4-911e-ff7fc254ba06.A07.gff')
# counts.dt <- counts.dt[counts.dt$V3 == 'exon', ]
# counts.dt$inc <- str_split_fixed(counts.dt$V9, '; ', 4)[, 3]
# counts.dt$inc <- gsub('inc "', '', counts.dt$inc)
# counts.dt$inc <- gsub('"', '', counts.dt$inc)
# counts.dt$inc <- as.numeric(counts.dt$inc)
# counts.dt$name <- c(1:nrow(counts.dt))
# exons <- counts.dt[, c('V1', 'V4', 'V5', 'name', 'inc', 'V7'), with = F]
# write.table(exons, '../exon_counts.bed', sep = '\t', row.names = F, col.names = F, quote = F)
files.K562 <- c('ENCFF636QII.A07.gff', 'ENCFF486POD.A07.gff')
files.HepG2 <- c('ENCFF670LIE.A07.gff', 'ENCFF074BOV.A07.gff')
path.to.folder <- '../exons/'
ParseFiles <- function(path.to.folder, file){
  dt <- fread(paste0(path.to.folder, file))
  dt <- subset(dt, V3 == 'exon')
  dt <- subset(dt, grepl('chr', dt$V1))
  dt$inc <- str_split_fixed(dt$V9, '; ', 4)[, 3]
  dt$inc <- gsub('"', '', dt$inc)
  dt$inc <- gsub('inc', '', dt$inc)
  dt$inc <- gsub(' ', '', dt$inc)
  dt$inc <- as.numeric(dt$inc)
  
  dt$exc <- str_split_fixed(dt$V9, '; ', 3)[, 2]
  dt$exc <- gsub('"', '', dt$exc)
  dt$exc <- gsub('exc', '', dt$exc)
  dt$exc <- gsub(' ', '', dt$exc)
  dt$exc <- as.numeric(dt$exc)
  
  dt <- dt[, c('V1', 'V4', 'V5', 'inc', 'exc', 'V7'), with = F]
  #colnames(dt) <- c('chr', 'start', 'stop', 'score', 'strand')
  dt
}

dt.list <- lapply(files.K562, function(file) ParseFiles(path.to.folder, file))
dt <- Reduce(function(...) merge(..., by = c('V1', 'V4', 'V5', 'V7'), all = T), dt.list)
dt[is.na(dt$inc.x), ]$inc.x <- 0 
dt[is.na(dt$inc.y), ]$inc.y <- 0 
dt[is.na(dt$exc.x), ]$exc.x <- 0 
dt[is.na(dt$exc.y), ]$exc.y <- 0 
dt$score <- dt$inc.x + dt$exc.y
dt$name <- c(1:nrow(dt))
dt.inc <- dt[, c('V1', 'V4', 'V5', 'name', 'score', 'V7'), with = F]
#dt.inc <- subset(dt.inc, score >= 10)
write.table(dt.inc, paste0(path.to.folder, 'exons_K562.bed'), sep = '\t', row.names = F, col.names = F, quote = F)

dt.psi <- dt
dt.psi$inc <- dt.psi$inc.x + dt.psi$inc.y
dt.psi$exc <- dt.psi$exc.x + dt.psi$exc.y
dt.psi <- subset(dt.psi, inc + 2*exc >= 20)
dt.psi$psi <- dt.psi$inc / (dt.psi$inc + 2 * dt.psi$exc)
dt.psi <- dt.psi[, c('V1', 'V4', 'V5', 'name', 'psi', 'V7'), with = F]
write.table(dt.psi, paste0(path.to.folder, 'exons_psi_K562.bed'), sep = '\t', row.names = F, col.names = F, quote = F)


# ends and starts of transcripts
path.to.anno <- '../conservative_features/gencode.v19.annotation.gtf'
dt <- fread(path.to.anno, skip = 5)
dt$transcript.id <- str_split_fixed(dt$V9, pattern = ';', Inf)[, 2] 
dt.CDS <- dt[dt$V3 == 'CDS', ]
CDS.transcript.id <- dt.CDS$transcript.id
dt <- dt[dt$transcript.id %in% CDS.transcript.id, ]
dt <- dt[dt$V3 == 'transcript', ]
dt$name <- c(1:nrow(dt))
dt$score <- 1
dt.ends <- dt[, c('V1', 'V5', 'V5', 'name', 'score', 'V7'), with = F]
dt.ends$end <- dt.ends$V5 + 1
dt.ends <- dt.ends[, c('V1', 'V5', 'end', 'name', 'score', 'V7'), with = F]
dt.starts <- dt[, c('V1', 'V4', 'V4', 'name', 'score', 'V7'), with = F]
dt.starts$end <- dt.starts$V4 + 1
dt.starts <- dt.starts[, c('V1', 'V4', 'end', 'name', 'score', 'V7'), with = F]
write.table(dt.ends, '../conservative_features/gencode.v19.annotation_transcripts_ends.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)
write.table(dt.starts, '../conservative_features/gencode.v19.annotation_transcripts_starts.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)

# exons no counts
dt <- fread(path.to.anno, skip = 5)
dt$exon.id <- str_split_fixed(dt$V9, pattern = ';', Inf)[, 10] 
dt.CDS <- dt[dt$V3 == 'CDS', ]
CDS.exon.id <- dt.CDS$exon.id
dt <- dt[dt$V3 == 'exon', ]
dt <- dt[dt$exon.id %in% CDS.exon.id, ]
dt$name <- c(1:nrow(dt))
dt$score <- 1
dt <- dt[, c('V1', 'V4', 'V5', 'name', 'score', 'V7'), with = F]
dt <- dt[!duplicated(dt[, c('V1', 'V4', 'V5', 'V7'), with = F]), ] 
write.table(dt, '../conservative_features/gencode.v19.annotation_exons_coding.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)

# poly A counts
fwd <- data.table(read.delim(pipe("grep \'^chr\' ../polyA/kidney_fwd"), header = F))
fwd <- fwd[c(grep('chr[1-9][1-2]?', fwd$V1), grep('chrX', fwd$V1), grep('chrY', fwd$V1)), ]
fwd$strand <- '+'
rev <- data.table(read.delim(pipe("grep \'^chr\' ../polyA/kidney_rev"), header = F))
rev <- rev[c(grep('chr[1-9][1-2]?', rev$V1), grep('chrX', rev$V1), grep('chrY', rev$V1)), ]
rev$strand <- '-'
res <- rbind(fwd, rev)
res$name <- c(1:nrow(res))
res <- res[, c('V1', 'V2', 'V3', 'name', 'V4', 'strand'), with = F]
write.table(res, '../polyA/kidney_all.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# poly A counts ENCODE
dt <- fread('../polyA/GSE30198_human.pas.uniq.bed')
dt <- subset(dt, grepl('AATAAA', dt$V4))
dt <- subset(dt, V5 > 50)
write.table(dt, '../polyA/GSE30198_human.pas.uniq_filtered.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# CAGE counts
fwd <- data.table(read.delim(pipe("grep \'^chr\' ../CAGE/NHEK_polyA+_fwd"), header = F))
fwd <- fwd[c(grep('chr[1-9][1-2]?', fwd$V1), grep('chrX', fwd$V1), grep('chrY', fwd$V1)), ]
fwd$strand <- '+'
rev <- data.table(read.delim(pipe("grep \'^chr\' ../CAGE/NHEK_polyA+_rev"), header = F))
rev <- rev[c(grep('chr[1-9][1-2]?', rev$V1), grep('chrX', rev$V1), grep('chrY', rev$V1)), ]
rev$strand <- '-'
res <- rbind(fwd, rev)
res$name <- c(1:nrow(res))
res <- res[, c('V1', 'V2', 'V3', 'name', 'V4', 'strand'), with = F]
write.table(res, '../CAGE/NHEK_polyA+_all.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# eCLIP
dt <- fread('../eClip/peaks_merged.bed')
dt$V4 <- paste(dt$V4, c(1:nrow(dt)), sep = '_')
dt$V5 <- 1
write.table(dt, '../eClip/peaks_merged_pretty.bed', sep = '\t', row.names = F, col.names = F, quote = F)


# evalues
evals <- fread('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/R_scape_estended_all.tsv')
evals <- evals[, c('id', 'E-value'), with = F]
evals$id <- str_split_fixed(evals$MSA, '_', 2)[, 2]
write.table(evals, '../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/R_scape_estended_all_pretty.tsv', 
            sep = '\t', row.names = F, col.names = T, quote = F)

# cryptic ss
files.K562 <- c('ENCFF636QII.A06.ssj.tsv', 'ENCFF486POD.A06.ssj.tsv')
files.HepG2 <- c('ENCFF670LIE.A06.ssj.tsv', 'ENCFF074BOV.A06.ssj.tsv')
path.to.folder <- '../cryptic_ss/'
ParseFiles <- function(path.to.folder, file){
  dt <- fread(paste0(path.to.folder, file))
  dt <- subset(dt, V6 == 'GTAG')
  dt$chr <- str_split_fixed(dt$V1, '_', 4)[, 1]
  dt$start <- as.numeric(str_split_fixed(dt$V1, '_', 4)[, 2])
  dt$stop <- as.numeric(str_split_fixed(dt$V1, '_', 4)[, 3])
  dt$strand <- str_split_fixed(dt$V1, '_', 4)[, 4]  
  dt.starts <- dt[, c('chr', 'start', 'strand', 'V2'), with = F]
  colnames(dt.starts) <- c('chr', 'coord', 'strand', 'V2')
  dt.stops <- dt[, c('chr', 'stop', 'strand', 'V2'), with = F]
  colnames(dt.stops) <- c('chr', 'coord', 'strand', 'V2')
  dt <- rbind(dt.starts, dt.stops)
  dt <- aggregate(V2 ~ chr + coord + strand, dt, sum)
  dt
}

dt.list <- lapply(files.K562, function(file) ParseFiles(path.to.folder, file))
dt <- Reduce(function(...) merge(..., by = c('chr', 'coord', 'strand'), all = T), dt.list)
dt[is.na(dt$V2.x), ]$V2.x <- 0 
dt[is.na(dt$V2.y), ]$V2.y <- 0 
dt$score <- dt$V2.x + dt$V2.y
dt <- dt[, c('chr', 'coord', 'strand', 'score')]
dt$expressed <- T

path.to.anno <- '../conservative_features/gencode.v19.annotation.gtf'
anno <- fread(path.to.anno, skip = 5)
anno <- subset(anno, V3 == 'exon')
anno.starts <- anno[, c('V1', 'V4', 'V7'), with = F]
colnames(anno.starts) <- c('chr', 'coord', 'strand')
anno.stops <- anno[, c('V1', 'V5', 'V7'), with = F]
colnames(anno.stops) <- c('chr', 'coord', 'strand')
anno <- rbind(anno.starts, anno.stops)
anno <- anno[!duplicated(anno), ]
anno$annotated <- T

ss <- merge(dt, anno, by = c('chr', 'coord', 'strand'), all = T)
ss[is.na(ss$expressed), ]$expressed <- F
ss[is.na(ss$annotated), ]$annotated <- F
cryptic.ss <- subset(ss, !annotated & expressed)
cryptic.ss$end <- cryptic.ss$coord + 1
cryptic.ss$name <- c(1:nrow(cryptic.ss))
cryptic.ss <- cryptic.ss[, c('chr', 'coord', 'end', 'name', 'score', 'strand')]
write.table(cryptic.ss, '../cryptic_ss/K562_cryptic_ss.bed', sep = '\t', row.names = F, col.names = T, quote = F)

# cryptic ss GTex
dt <- fread('../cryptic_ss/cryptic_sites_GTEx.tsv')
dt$end <- dt$pos + 1
dt$counts <- 1
dt$id <- c(1:nrow(dt))
dt <- dt[, c('chr', 'pos', 'end', 'id', 'counts', 'strand'), with = F]
write.table(dt, '../cryptic_ss/cryptic_sites_GTEx.bed', sep = '\t', row.names = F, col.names = F, quote = F)


# circle RNA CIRCpedia
dt <- fread('../circRNA/CIRCpedia.csv')
colnames(dt) <- c('ID', 'Species', 'Gene', 'Transcript', 'Location', 'strand', 
                  'FPM', 'ExonStart_ExonEnd', 'Seq_type', 'Cell_line', 'Conservation', 
                  'MapsSplice', 'Enrichment')
write.table(dt, '../circRNA/CIRCpedia.csv', sep = '\t', row.names = F, quote = F)
dt <- subset(dt, Cell_line == 'K562')
dt$chr <- str_split_fixed(dt$Location, ':', 2)[, 1]
dt$Location <- str_split_fixed(dt$Location, ':', 2)[, 2]
dt$start <- as.numeric(str_split_fixed(dt$Location, '-', 2)[, 1])
dt$stop <- as.numeric(str_split_fixed(dt$Location, '-', 2)[, 2])
dt <- dt[, c('chr', 'start', 'stop', 'ID', 'FPM', 'strand'), with = F]
write.table(dt, '../circRNA/CIRCpedia_K562.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# circle RNA TSCD
dt <- fread('../circRNA/TSCD.txt')
colnames(dt) <- c('Tissue', 'Coordinates', 'chr', 'Donor_site', 'Acceptor_site', 'Junction_reads', 'strand',
                  'Algorithms', 'SRPTM', 'Gene_annotation', 'Gene_body', 'Gene_type', 'gene_strand', 'microRNAs', 
                  'RBPs', 'Flanking_length(donor to acceptor site)')
write.table(dt, '../circRNA/TSCD.txt', sep = '\t', row.names = F, quote = F)
dt <- subset(dt, grepl('liver', Tissue))
dt <- dt[, c('Tissue', 'chr', 'Donor_site', 'Acceptor_site', 'strand', 'SRPTM', 'Coordinates'), with = F]
dt$counts <- unlist(lapply(dt$SRPTM, function(value) mean(as.numeric(unlist(strsplit(value, ','))))))
tmp <- subset(dt, grepl(',', strand))
dt <- subset(dt, !grepl(',', strand))
tmp$strand <- '-'
dt <- rbind(dt, tmp)
tmp$strand <- '+'
dt <- rbind(dt, tmp)
dt$name <- c(1:nrow(dt))
dt <- dt[, c('chr', 'Donor_site', 'Acceptor_site', 'name', 'counts', 'strand'), with = F]
write.table(dt, '../circRNA/TSCD_liver.bed', sep = '\t', row.names = F, col.names = F, quote = F)
