library(data.table)
library(ggplot2)
options(scipen = 999)
colours <- c('0,0,0', '0,0,255', '255,0,0')
path <- 'data/out/'
ph <- fread('../conservative_features/whole_human_genome/alignments_whole_human_genome_processed_norm_structures.tsv')
oncogenes <- fread('../ongene_human.txt')
flank_length <- 1000

FilterAndMakeTrack <- function(tissue, feature, path, ph, oncogenes, flank_length){
  dt <- fread(paste0(path, tissue, '_', feature, 'filtered_denominator30_qval_qsign.tsv'))
  dt <- dt[dt$q < 0.05, ]
  if(feature == 'site'){
    dt$start <- dt$coord
    dt$end <- dt$coord + 1    
  }
  dt$name <- paste0('deltaPSI=', dt$delta_psi_c, "|inc_normal=", dt$inc_normal, 
                    "|ret_normal=", dt$exc_normal, "|inc_tumor=", dt$inc_tumor, 
                    "|ret_tumor=", dt$exc_tumor)
  dt$score <- '1'
  dt$itemRGB <- ''
  dt[dt$delta_psi_c < 0,]$itemRGB <- colours[2]
  dt[dt$delta_psi_c >= 0,]$itemRGB <- colours[3]
  dt$start <- dt$start
  dt$thickStart <- dt$start
  dt$thickEnd <- dt$end
  bed <- dt[, c('chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRGB')]
  
  ph <- ph[, c('id', 'panhandle_start', 'panhandle_end', 'panhandle_left_hand', 'panhandle_right_hand', 'chr'), with = F]
  dt <- dt[dt$gene.name %in% oncogenes$OncogeneName, ]
  merged <- merge(dt, ph, by = 'chr', all = F, allow.cartesian = TRUE)
  merged <- merged[abs(merged$panhandle_start - merged$start) <= flank_length | abs(merged$panhandle_left_hand - merged$start) <= flank_length |
                     abs(merged$panhandle_right_hand - merged$start) <= flank_length | abs(merged$panhandle_end - merged$start) <= flank_length]
  merged$coords <- paste0(merged$chr, ":", merged$start, "-", merged$end)
  merged <- merged[! duplicated(merged$start), ]
  merged <- merged[, names(merged) %in% c('coords', 'name', 'gene.name'), with = F]
  
  print(paste0(tissue, " ", feature))
  print(merged)
  write.table(merged, paste0(path, "/", tissue, "_", feature, "_", "with_ph.tsv"), sep = "\t", quote = F, row.names = F)
  track_name <- paste0(tissue, '_', feature, '_', 'oncogenes')
  file_name <- paste0(path, tissue, '_', feature, '_track.bed') 
  f <- file(file_name, "w")
  writeLines(paste0('track name="', track_name, '"', ' itemRgb="On"'),f)
  write.table(bed, f, sep = '\t', row.names = F, col.names = F, quote = F)  
  close(f)
  system(paste0('bedtools sort -i ', file_name, ' > track_sorted.bed'))
  system(paste0('../tools/bedToBigBed track_sorted.bed ../tools/hg19.chrom.sizes ', gsub('.bed', '.bb', file_name)))
}

FilterAndMakeTrack('kidney', 'site', path, ph, oncogenes, flank_length)
FilterAndMakeTrack('liver', 'site', path, ph, oncogenes, flank_length)
FilterAndMakeTrack('kidney', 'intron', path, ph, oncogenes, flank_length)
FilterAndMakeTrack('liver', 'intron', path, ph, oncogenes, flank_length)

###############################
dt <- fread('data/kidney_site.tsv')
dt.tmp <- dt[dt$chr == 'chr1' & dt$coord == 110884890, grep('D', colnames(dt), value = T), with = F]
dt.tmp[is.na(dt.tmp)] <- 0
donors <- unique(gsub("counts_ret_tumor", "", grep("counts_ret_tumor", colnames(dt.tmp), value = T)))
for(donor in donors){
  set(dt.tmp, j  = paste0("denom_normal_", donor), value = dt.tmp[,paste0("counts_ret_normal", donor),with = F] + dt.tmp[,paste0("counts_inc_normal", donor),with = F]) 
  set(dt.tmp, j  = paste0("denom_tumor_", donor), value =  dt.tmp[,paste0("counts_ret_tumor", donor),with = F] + dt.tmp[,paste0("counts_inc_tumor", donor),with = F]) 
  
  set(dt.tmp, j  = paste0("psi_normal_", donor), value = dt.tmp[,paste0("counts_inc_normal", donor),with = F] / (dt.tmp[,paste0("counts_ret_normal", donor),with = F] + dt.tmp[,paste0("counts_inc_normal", donor),with = F])) 
  set(dt.tmp, j  = paste0("psi_tumor_", donor), value =  dt.tmp[,paste0("counts_inc_tumor", donor),with = F] / (dt.tmp[,paste0("counts_ret_tumor", donor),with = F] + dt.tmp[,paste0("counts_inc_tumor", donor),with = F])) 
  set(dt.tmp, j  = paste0("delta_psi_", donor), value =  dt.tmp[,paste0("psi_tumor_", donor),with = F] - dt.tmp[,paste0("psi_normal_", donor),with = F])
  
}
dt.tmp <- dt.tmp[, grep('psi', colnames(dt.tmp), value = T), with = F]# denom # psi # delta psi
to.plot <- as.data.frame(t(dt.tmp))
to.plot$type <- 'tumor'
to.plot[grep('normal', row.names(to.plot)), ]$type <- "normal"
to.plot <- to.plot[!is.na(to.plot$V1), ]
png('data/kidney_site_RBM15_psi.png', width = 6, height = 6, units = 'in', res = 300)
ggplot(to.plot, aes(x = type, y = V1, fill = type)) +
  geom_boxplot() +
  xlab(label = "") +
  ylab(label = "psi") +
  ggtitle("RBM15")
dev.off()
