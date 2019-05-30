options(scipen = 999)
library(data.table)
library(dplyr)
library(stringr)
MakeFlankedHandles <- function(ph.handle.path.bed, flank.length = 1000){
  handles <- fread(ph.handle.path.bed)
  handles$V2 <- handles$V2 - flank.length
  handles$V3 <- handles$V3 + flank.length
  write.table(handles, gsub('.bed', '_with_flanks.bed', ph.handle.path.bed), row.names = F, col.names = F, quote = F, sep = '\t')
}

IntersectPHandMut <- function(path.to.mut, ph.handle.path.bed, tissue, dataset){
  x <- paste0('awk -F"\\t" \'{print $9,$10,$11,$1"_"$2"_"$15"_"$16"_"$17}\' OFS="\\t" ', path.to.mut, ' | sed \'s/^/chr/\' | tail -n +2 |sort -u | bedtools sort -i stdin > tmp_mut.bed')
  system(x)
  path.to.mut <- 'tmp_mut.bed'
  x <- paste0('bedtools sort -i ', ph.handle.path.bed, ' | bedtools intersect -a stdin -b ', path.to.mut, ' -wa -wb')
  mut.in.df <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  mut.in.df$ph.id <- as.numeric(str_split_fixed(mut.in.df$V4, '_', 2)[, 1])
  mut.in.df$handle.w.mut <- str_split_fixed(mut.in.df$V4, '_', 2)[, 2]
  mut.in.df <- mut.in.df[, c('ph.id', 'handle.w.mut', 'V8')]
  colnames(mut.in.df) <- c('ph.id', 'handle.w.mut', 'mut.id')
  write.table(mut.in.df, paste0('data/', tissue, '_', dataset, '_ph_with_mut.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
}

IntersectPanhandles <- function(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out, ph.with.mut){
  feature.path <- paste0(prefix, tissue, '_', feature, '_filtered_population.tsv')
  print(feature.path)
  feature.dt <- fread(feature.path)
  if(feature == 'exon' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)    
  } else if(feature == 'exon' & looped.out == F){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)   
  } else if(grepl('site', feature)){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x)   
  }
  else if(feature == 'intron' & looped.out == F){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x) 
    dt <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
    ncol <- dim(feature.dt)[2]
    dt <- dt[dt$V2 <= dt[[ncol + 3]] & dt$V2 >= dt[[ncol + 2]] | dt$V3 <= dt[[ncol + 3]] & dt$V3 >= dt[[ncol + 2]], ]
    write.table(dt, paste0(prefix, tissue, '_', feature, '_overlapped.tsv'), row.names = F, col.names = F, quote = F, sep = '\t')
  } else if(feature == 'intron' & looped.out){
    x <- paste0('tail -n +2 ', feature.path, ' | bedtools sort -i stdin', ' | bedtools intersect -a stdin -b ', ph.handle.flanked.path.bed, ' -wa -wb > ', prefix, tissue, '_', feature, '_overlapped.tsv')
    system(x) 
    dt <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
    ncol <- dim(feature.dt)[2]
    dt <- dt[dt$V2 <= dt[[ncol + 3]] & dt$V2 >= dt[[ncol + 2]] | dt$V3 <= dt[[ncol + 3]] & dt$V3 >= dt[[ncol + 2]], ]
    dt$ph.id <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 1]
    dt$handle <- str_split_fixed(dt[[ncol + 4]], '_', 2)[, 2]
    overlapped <- merge(dt, dt, by = c(paste("V", c(1:ncol), sep = ''), 'ph.id'))
    overlapped <- overlapped[overlapped$handle.y != overlapped$handle.x, ]
    overlapped$V21.x <- overlapped$V21.x + 1000
    overlapped$V21.y <- overlapped$V21.y + 1000
    overlapped$V22.x <- overlapped$V22.x - 1000
    overlapped$V22.y <- overlapped$V22.y - 1000
    overlapped <- overlapped[overlapped$V2 <= overlapped$V21.x & 
                               overlapped$V2 <= overlapped$V21.y & 
                               overlapped$V3 >= overlapped$V22.x & 
                               overlapped$V3 >= overlapped$V22.y,]
    to.drop <- grep('x|y', colnames(overlapped), value = T)
    overlapped <- overlapped[, ! names(overlapped) %in% c(to.drop), with = F]
    write.table(overlapped, paste0(prefix, tissue, '_', feature, '_overlapped.tsv'), row.names = F, col.names = T, quote = F, sep = '\t')
  }
  overlapped <- fread(paste0(prefix, tissue, '_', feature, '_overlapped.tsv'))
  ncol <- dim(feature.dt)[2]
  overlapped$overlap <- T
  if(!any(grepl('ph.id', colnames(overlapped)))){
    overlapped$ph.id <- as.numeric(str_split_fixed(overlapped[[ncol + 4]], '_', 2)[, 1])  
  }
  dt <- merge(feature.dt, overlapped[, c('V1', 'V2', 'V3', 'V4', 'overlap', 'ph.id'), with = F], 
              by.x = c('chr', 'start', 'end', 'strand'), 
              by.y = c('V1', 'V2', 'V3', 'V4'), all.x = T)
  dt[is.na(dt$overlap), ]$overlap <- F
  mut.dt <- fread(ph.with.mut)
  dt <- merge(dt, mut.dt, by = 'ph.id', all.x = T)
  dt$ph.w.mut <- F
  dt[!is.na(dt$mut.id), ]$ph.w.mut <- T
  if(grepl('site', feature)){
    dt <- dt[dt$end - dt$start <= 500, ]
  }
  dt$feature.id <- paste(dt$chr, dt$start, dt$end, dt$strand, sep = '_')
  dt
}



TestMut <- function(tissue, feature, locations, ph.path.bed, ph.handle.flanked.path.bed, prefix, ph.with.mut, dataset, programme){
  results <- list()
  counts <- list()
  dt.nlooped.out <- IntersectPanhandles(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out = F, ph.with.mut)
  write.table(dt.nlooped.out, paste0('data/', programme, '_', dataset, '_', tissue, '_', feature, '_close_with_all_info.tsv'), sep = '\t', row.names = F, quote = F)
  counts[['total_denom_filtered']] <- length(unique(dt.nlooped.out$feature.id))
  counts[['DE']] <- length(unique(dt.nlooped.out[dt.nlooped.out$significant, ]$feature.id))
  counts[['close_struct']] <- length(unique(dt.nlooped.out[dt.nlooped.out$overlap, ]$feature.id))
  counts[['close_struct_and_DE']] <- length(unique(dt.nlooped.out[dt.nlooped.out$significant & dt.nlooped.out$overlap, ]$feature.id))
  counts[['close_struct_with_mut']] <- length(unique(dt.nlooped.out[dt.nlooped.out$overlap & dt.nlooped.out$ph.w.mut, ]$feature.id))
  counts[['close_struct_with_mut_and_DE']] <- length(unique(dt.nlooped.out[dt.nlooped.out$significant & dt.nlooped.out$overlap & dt.nlooped.out$ph.w.mut, ]$feature.id))
  interesting.cases.nlooped.out <- dt.nlooped.out[dt.nlooped.out$significant & dt.nlooped.out$overlap & dt.nlooped.out$ph.w.mut, c('chr', 'start', 'end', 'ph.id'), with = F]
  interesting.cases.nlooped.out <- interesting.cases.nlooped.out[!duplicated(interesting.cases.nlooped.out[, c('chr', 'start', 'end'), with = F]), ]
  interesting.cases.looped.out <- data.table(chr = character(), start = numeric(), end = numeric(), ph.id = numeric())
  counts[['looped.out']] <- 0
  counts[['looped.out_and_DE']] <- 0
  counts[['looped.out_with_mut']] <- 0
  counts[['looped.out_with_mut_and_DE']] <- 0
  if(length(locations) > 1){
    dt.looped.out <- IntersectPanhandles(ph.path.bed, ph.handle.flanked.path.bed, prefix, feature, tissue, looped.out = T, ph.with.mut)
    write.table(dt.nlooped.out, paste0('data/', programme, '_', dataset, '_', tissue, '_', feature, '_looped_out_with_all_info.tsv'), sep = '\t', row.names = F, quote = F)
    counts[['looped.out']] <- length(unique(dt.looped.out[dt.looped.out$overlap, ]$feature.id))
    counts[['looped.out_and_DE']] <- length(unique(dt.looped.out[dt.looped.out$significant & dt.looped.out$overlap, ]$feature.id))
    counts[['looped.out_with_mut']] <- length(unique(dt.looped.out[dt.looped.out$overlap & dt.looped.out$ph.w.mut, ]$feature.id))
    counts[['looped.out_with_mut_and_DE']] <- length(unique(dt.looped.out[dt.looped.out$significant & dt.looped.out$overlap & dt.looped.out$ph.w.mut, ]$feature.id))
    interesting.cases.looped.out <- dt.looped.out[dt.looped.out$significant & dt.looped.out$overlap & dt.looped.out$ph.w.mut, c('chr', 'start', 'end', 'ph.id'), with = F]
    interesting.cases.looped.out <- interesting.cases.looped.out[!duplicated(interesting.cases.looped.out[, c('chr', 'start', 'end'), with = F]), ]
  }
  counts.df <- as.data.frame(counts)
  counts.df$tissue <- tissue
  counts.df$feature <- feature
  counts.df$dataset <- dataset
  counts.df$programme <- programme
  
  interesting.cases.looped.out$location <- 'looped.out'
  interesting.cases.nlooped.out$location <- 'close'
  
  interesting.cases <- rbind(interesting.cases.looped.out, interesting.cases.nlooped.out)
  interesting.cases$tissue <- tissue
  interesting.cases$feature <- feature
  interesting.cases$dataset <- dataset
  interesting.cases$programme <- programme
  
  results[['interesting.cases']] <- interesting.cases
  results[['counts']] <- counts.df
  results
}

MakeCounts <- function(x){
  programmes <- c('spladder', 'ipsa')
  datasets <- c('whole', 'filtered')
  tissues <- c('liver', 'kidney')
  counts.list <- list()
  interesting.cases <- list()
  for(programme in programmes){
    print(programme)
    if(programme == 'spladder'){
      features <- c('exon', 'intron', '5prime_site', '3prime_site')
      out.path <- 'data/spladder_output/'
      prefix <- 'data/spladder_output/spladder_'
    } else if(programme == 'ipsa'){
      features <- c('exon', 'intron', 'site')
      out.path <- 'data/ipsa_output/'
      prefix <- 'data/ipsa_output/ipsa_'
    }
    for(dataset in datasets){
      print(dataset)
      if(dataset == 'filtered'){
        ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered_handles.bed'
        ph.path.bed <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_filtered.bed'
        out.path2 <- paste0(out.path, 'filtered_ds_pictures/')
      } else if(dataset == 'whole'){
        ph.handle.path.bed <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed_handles.bed'
        ph.path.bed <- '../python_scripts/folding_pretty_copy/out/panhandles_preprocessed.bed'
        out.path2 <- paste0(out.path, 'whole_ds_pictures/')
      }
      ph.handle.flanked.path.bed <- gsub('.bed', '_with_flanks.bed', ph.handle.path.bed)
      if (!file.exists(ph.handle.flanked.path.bed)){
        print('Making flanked handles..')
        MakeFlankedHandles(ph.handle.path.bed)
      }
      print('I have flanked handles file!')
      for(tissue in tissues){
        print(tissue)
        if(tissue == 'kidney'){
          path.to.mut <- '../mutations_icgc/SNVs/RECA_simple_somatic_mutation.open_only_SNV.tsv'
        } else if(tissue == 'liver'){
          path.to.mut <- '../mutations_icgc/SNVs/LIRI_simple_somatic_mutation.open_only_SNV_hypermut_donot_removed.tsv'
        }
        ph.with.mut <- paste0('data/', tissue, '_', dataset, '_ph_with_mut.tsv')
        if(!file.exists(ph.with.mut)){
          print('Making mutation file..')
          IntersectPHandMut(path.to.mut, ph.handle.path.bed, tissue, dataset)
        }
        print('I have mutation file!')
        for(feature in features){
          print(feature)
          if(grepl('site', feature)){
            locations <- c(FALSE)
          } else{
            locations <- c(FALSE, TRUE)
          }
          experiment <- paste0(programme, '_', dataset, '_', tissue, '_', feature)
          print(experiment)
          results <- TestMut(tissue, feature, locations, ph.path.bed, ph.handle.flanked.path.bed, prefix, ph.with.mut, dataset, programme)
          counts.list[[experiment]] <- results[['counts']]
          interesting.cases[[experiment]] <- results[['interesting.cases']]
        }
      }
    }
  }
  interesting.cases.df <- do.call("rbind", interesting.cases) 
  interesting.cases.df <- interesting.cases.df[!duplicated(interesting.cases.df[, c('chr', 'start', 'end'), with = F]), ]
  for(exp in names(counts.list)){
    dt <- counts.list[[exp]]
    if(!any(grepl('loop', colnames(dt)))){
      dt$looped.out <- 0
      dt$looped.out_and_DE <- 0
      dt$looped.out_with_mut <- 0
      dt$looped.out_with_mut_and_DE <- 0
      counts.list[[exp]] <- dt
    }
  }
  counts.df <- do.call('rbind', counts.list)
  write.table(counts.df, 'data/counts_all_features_ph_mut.tsv', sep = '\t', row.names = F, quote = F)
  write.table(interesting.cases.df, 'data/cases_all_features_ph_mut_2.tsv', sep = '\t', row.names = F, quote = F)
}

MakeCounts(1)

MatchMutAndFeature <- function(dt){
  dt <- fread('data/cases_most_interesting.csv', sep = ',')
  for(tissue. in unique(dt$tissue)){
    ph_mut <- fread(paste0('data/', tissue., '_filtered_ph_with_mut.tsv'))
    dt <- merge(dt, ph_mut[, c('ph.id', 'mut.id'), with = F], by = 'ph.id', all.x = T)
    if(tissue. == 'liver'){
      mut <- fread('../mutations_icgc/SNVs/LIRI_simple_somatic_mutation.open_only_SNV_hypermut_donot_removed.tsv')
    } else if(tissue. == 'kidney'){
      mut <- fread('../mutations_icgc/SNVs/RECA_simple_somatic_mutation.open_only_SNV.tsv')
    }
    dt$icgc_mutation_id <- str_split_fixed(dt$mut.id, '_', 2)[, 1]
    dt <- merge(dt, mut[, c('icgc_mutation_id', 'icgc_donor_id'), with = F], by = 'icgc_mutation_id', all.x = T)
    dt$DE <- F
    dt$exon_chr_start_end <- paste(dt$chr, dt$start, dt$end, sep = '_')
    dt <- dt[!duplicated(dt), ]
    colnames(dt)[colnames(dt) == 'icgc_donor_id'] <- 'icgc_donor_id_mut'
    for(programme. in unique(dt$programme)){
      for(feature. in unique(dt$feature)){
        if(programme. == 'spladder'){
          features.dt <- fread(paste0('data/spladder_output/merge_graphs_', tissue., '_', feature., '.txt.gz'))
        } else if(programme. == 'ipsa'){
          features.dt <- fread(paste0('data/ipsa_output/output_', tissue., '_', feature., '.tsv'))
          if(!grepl('chr', colnames(features))){
            features.dt$chr <- str_split_fixed(features.dt$exon_chr_start_end_strand, ':', 2)[, 1]
            features.dt$coord <- str_split_fixed(features.dt$exon_chr_start_end_strand, ':', 2)[, 2]
            features.dt$start <- as.numeric(str_split_fixed(features.dt$coord, '_', 2)[, 1])
            features.dt$end <- str_split_fixed(features.dt$coord, '_', 2)[, 2]
            features.dt$end <- as.numeric(unlist(lapply(features.dt$end, function(x){substr(x, 1, nchar(x) - 1)})))
          } else if(grepl('coord', colnames(features.dt))){
            features.dt$start <- features.dt$coord
            features.dt$end <- features.dt$coord + 1
          }
          features.dt$exon_chr_start_end <- paste(features.dt$chr, features.dt$start, features.dt$end, sep = '_')
          features.dt <- features.dt[features.dt$exon_chr_start_end %in% 
                                       dt[dt$programme == programme. & dt$tissue == tissue. & dt$feature == feature. ,]$exon_chr_start_end, ]
          for(exon_chr_start_end. in unique(dt[dt$programme == programme. & dt$tissue == tissue. & dt$feature == feature. ,]$exon_chr_start_end)){
            for(donor. in unique(dt[dt$programme == programme. & dt$tissue == tissue. & dt$feature == feature. & dt$exon_chr_start_end == exon_chr_start_end. ,]$icgc_donor_id_mut)){
              features.tmp <- features.dt[features.dt$exon_chr_start_end == exon_chr_start_end., grep(donor., colnames(features.dt), value = T), with = F]
              if(dim(features.tmp)[1] == 0){
                dt[dt$programme == programme. & dt$tissue == tissue. & dt$feature == feature. & dt$exon_chr_start_end == exon_chr_start_end. & dt$icgc_donor_id_mut == donor., ]$DE <- F
              } else{
                 
              }
            }
          }
          
        } 
      }
      
    }
  }
  
  
}




