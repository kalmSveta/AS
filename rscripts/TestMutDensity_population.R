#!/usr/bin/Rscript
.libPaths(c( .libPaths(), "../R/x86_64-redhat-linux-gnu-library/3.4/"))
library(stringr)
library(stringi)
library(data.table)
options(scipen = 999)

# TestMutDensity.R path.to.sample path.to.population path.to.mut.folder.sample path.to.mut.folder.population path.to.phs count.patients count.stacked subtract
args <- commandArgs(trailingOnly = TRUE)
path.to.sample <- args[1]
path.to.population <- args[2]
path.to.mut.folder.sample <- args[3]
path.to.mut.folder.population <- args[4]
path.to.phs <- args[5]
count.patients <- as.logical(args[6])
count.stacked <- as.logical(args[7])
subtract <- as.logical(args[8])




CountLength <- function(path.to.gr, add){
  x <- paste0('bedtools sort -i ', path.to.gr, ' | bedtools merge -i stdin -s ', add, ' | awk -F"\\t" \'{print $1,$3-$2+1}\' OFS="\\t" | awk \'{s+=$2}END{print s}\'')
  length_ <- tryCatch(scan(pipe(x)), error = function(e) 0)
  length_ <- as.numeric(length_)
  print(length_)
  return(length_)
}

SelectStacked <- function(dt, path.to.phs){
  print('Select stacked..')
  phs <- fread(path.to.phs)
  phs$id <- as.character(phs$id)
  dt$id <- str_split_fixed(dt$V4, '_', 2)[, 1]
  dt$handle <- str_split_fixed(dt$V4, '_', 2)[, 2]
  dt <- merge(dt, phs[, c('id', 'structure', 'panhandle_start', 'panhandle_right_hand', 'alignment1', 'alignment2'), with = F], all.x = T)
  dt$local.mut.coord <- dt$V8 - dt$V2 + 1
  dt$handle.struct <- ''
  dt[dt$handle == 'left']$handle.struct <- unlist(lapply(dt[dt$handle == 'left']$structure, function(struct){
    pos <- stri_locate_all(pattern = ')', struct, fixed = T)[[1]][1, 'start']
    substr(struct, 1, (pos - 1))
  }))
  dt[dt$handle == 'right']$handle.struct <- unlist(lapply(dt[dt$handle == 'right']$structure, function(struct){
    pos <- stri_locate_all(pattern = ')', struct, fixed = T)[[1]][1, 'start']
    substr(struct, pos, nchar(struct))
  }))
  dt$mut.nt.struct <- unlist(apply(dt, 1, function(row){
    substr(row['handle.struct'], as.numeric(row['local.mut.coord']), as.numeric(row['local.mut.coord']))
  }))
  dt <- subset(dt, mut.nt.struct != '.')
  dt
}

CountMutations <- function(path.to.mut.folder, count.patients, count.stacked, path.to.phs){
  list.of.files <- file.info(list.files(path.to.mut.folder, full.names = T))
  list.of.non.empty.files <- subset(list.of.files, size != 0)
  counts <- unlist(lapply(rownames(list.of.non.empty.files), function(file.) {
    print(file.)
    dt <- fread(file.)
    # Count simple muattions or mutations in patients
    if(count.patients){
      print('Counting patients..')
      patient.counts <- unlist(apply(dt, 1, function(row) length(grep('1', row[16:ncol(dt)]))))
      dt$patient.counts <- patient.counts
    } else{
      dt$patient.counts <- 1
    }
    dt <- dt[, c(colnames(dt)[c(1:11)], 'patient.counts'), with = F]
    dt <- dt[!duplicated(dt[, c(7,8,9), with = F]), ]
    # Count mutations in handles or in stacked pairs
    if(count.stacked){
      dt <- SelectStacked(dt, path.to.phs)
    }
    sum(dt$patient.counts)
  }))
  sum(counts)
}

TestMutEnrichment <- function(path.to.sample, path.to.population, path.to.mut.folder.sample, path.to.mut.folder.population, path.to.phs,
                              count.patients, count.stacked, subtract = 'False'){
  # Count length of merged handles
  print('Count length sample..')
  add <- ""
  out <- CountLength(path.to.sample, add)
  if(count.patients){
    length1 <- out * 2504  
  } else{
    length1 <- out
  }
  
  
  # Count length of merged conin
  print('Count length population..')
  if(subtract == 'True'){
    print('Subtracting sample from popuation..')
    add <- paste0("| bedtools subtract -s -a stdin -b ", path.to.sample)
  }
  out <- CountLength(path.to.population, add)
  if(count.patients){
    length2 <- out * 2504  
  } else{
    length2 <- out
  }
  
  # Count mutations
  count.mut1 <- CountMutations(path.to.mut.folder.sample, count.patients, count.stacked, path.to.phs)
  count.mut2 <- CountMutations(path.to.mut.folder.population, count.patients, FALSE, path.to.phs)
  
  
  print(paste0('sample: length=', as.character(length1), ', mutations=', as.character(count.mut1)))
  print(paste0('population: length=', as.character(length2), ', mutations=', as.character(count.mut2)))
  p.value <- poisson.test(x = c(count.mut1, count.mut2), T = c(length1, length2), alternative = 't', r = 1)$p.value
  return(list('proc1' = count.mut1 / length1,
              'proc2' = count.mut2 / length2,
              'p.value' = p.value))
}

print(TestMutEnrichment(path.to.sample, path.to.population, path.to.mut.folder.sample, path.to.mut.folder.population, path.to.phs,
                        count.patients, count.stacked, subtract))

