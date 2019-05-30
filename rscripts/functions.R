options(scipen = 999)
library(data.table)

Make_bed = function(df, file_name, different_colours = F, track_name){
  colours <- c('0,0,0','0,0,255', '255,0,0')
  #df <- df[order(df[, order_by]), ]
  #df$index <- paste(df$id, df$energy, sep='_')
  df$score <- df$energy / 100
  if(different_colours){
    names(colours) <- c(T, F)
    df$itemRGB <- colours[1]
    df[df$delta_psi_c < 0,]$itemRGB <- colours[2]
    df[df$delta_psi_c >= 0,]$itemRGB <- colours[3]
    add <- 'itemRgb="On"'
  } else {
    df$itemRGB <- '0,0,0'
    add <- ''
  }
  df$block_count <- 2
  df$blockSizes <- paste(as.character(df$al1_length),as.character(df$al2_length-1),sep=',')
  df$blockStarts <- paste(as.character(0),as.character(df$panhandle_right_hand - df$panhandle_start),sep=',')
  print(head(df))
  bed <- df[, c('chr','panhandle_start','panhandle_end','index','score','strand','panhandle_start','panhandle_end','itemRGB',
                'block_count','blockSizes','blockStarts')]
  colnames(bed) <- c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes',
                     'blockStarts')
  print(head(bed))
  f <- file(file_name, "w")
  writeLines(paste('track name="',track_name, '" ', add, sep=''),f)
  write.table(bed,f,sep = '\t',row.names = F,col.names = F,quote = F)  
  close(f)
  return(bed)
}

BedToBigBed <- function(path.to.ph.bed, n_header_lines = 0){
  x <- paste0('tail -n +', n_header_lines + 1, ' ', path.to.ph.bed, 
              ' | bedtools sort -i stdin > sorted.bed')
  system(x)
  x <- paste0('../tools/bedToBigBed sorted.bed ', '../tools/hg19.chrom.sizes ', gsub('bed', 'bb', path.to.ph.bed))
  system(x)
  system('rm sorted.bed')
}

MakeHandlesBed <- function(dt){
  dt.left <- dt[, c('chr', 'panhandle_start', 'panhandle_left_hand', 'id')]
  dt.right <- dt[, c('chr', 'panhandle_right_hand', 'panhandle_end', 'id')]
  colnames(dt.left) <- c('chr', 'start', 'stop', 'name')
  colnames(dt.right) <- c('chr', 'start', 'stop', 'name')
  dt.left$name <- paste(dt.left$name, '_left', sep = '')
  dt.right$name <- paste(dt.right$name, '_right', sep = '')
  dt.both <- rbind(dt.left, dt.right)
  dt.both
}

PreprocessAfterRecalcForRandom <- function(input.path, before.mut.path, panh.path){
  mutated <- fread(input.path, header = T)
  before.mut <- fread(before.mut.path, header = T)
  df <- fread(panh.path, header = T)
  
  mutated <- mutated[mutated$new_energy != 'new_energy',]
  mutated$new_end_al1 <- as.numeric(mutated$new_end_al1)
  mutated$new_end_al2 <- as.numeric(mutated$new_end_al2)
  mutated$new_start_al2 <- as.numeric(mutated$new_start_al2)
  mutated$new_start_al1 <- as.numeric(mutated$new_start_al1)
  mutated$new_al2_length <- mutated$new_end_al2 - mutated$new_start_al2 + 1
  mutated$new_al1_length <- mutated$new_end_al1 - mutated$new_start_al1 + 1
  mutated$new_energy <- as.numeric(mutated$new_energy)
  mutated$panhandle_id <- as.numeric(mutated$panhandle_id)
  
  zero.structures <- setdiff(unique(before.mut$panhandle_id), unique(mutated$panhandle_id))                                                                                                                        
  zero.ztructures.dt <- before.mut[before.mut$panhandle_id %in% zero.structures, c('panhandle_id', 'mut_coord')]
  zero.ztructures.dt$V1 <- 0
  zero.ztructures.dt$new_alignment1 <- ''
  zero.ztructures.dt$new_alignment2 <- ''                                                                                                                                                                   
  zero.ztructures.dt$new_end_al1 <- 0
  zero.ztructures.dt$new_end_al2 <- 0                                                                                                                                                                       
  zero.ztructures.dt$new_start_al2 <- 0                                                                                                                                                                     
  zero.ztructures.dt$new_start_al1 <- 0                                                                                                                                                                     
  zero.ztructures.dt$new_energy <- 0
  zero.ztructures.dt$new_al2_length <- 0
  zero.ztructures.dt$new_al1_length <- 0
  zero.ztructures.dt$new_energy <- 0
  print(dim(zero.ztructures.dt))
  mutated <- rbind(mutated, zero.ztructures.dt)
  mutated.merged <- merge(mutated, df, by.x = 'panhandle_id', by.y = 'id' )
  mutated.merged$dE <- mutated.merged$new_energy - mutated.merged$energy
  mutated.merged
}

Zscore <- function(a, A, b, B){
  pa <- a / A
  pb <- b / B
  z <- (pa - pb) / sqrt(pa * (1 - pa) / A + pb * (1 - pb) / B)
  z
}
