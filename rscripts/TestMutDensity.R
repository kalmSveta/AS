#!/usr/bin/Rscript
library(stringr)
options(scipen = 999)

# TestMutDensity.R path.to.sample path.to.population path.to.mut subtract
args <- commandArgs(trailingOnly = TRUE)
path.to.sample <- args[1]
path.to.population <- args[2]
path.to.mut <- args[3]
subtract <- args[4]


CountMutations <- function(path.to.gr, path.to.mut, add){
  x <- paste0('bedtools sort -i ', path.to.gr, ' | bedtools merge -i stdin ', add, ' | awk -F"\\t" \'{print $1,$3-$2+1}\' OFS="\\t" | awk \'{s+=$2}END{print s}\'')
  length_ <- tryCatch(scan(pipe(x)), error = function(e) 0)
  length_ <- as.numeric(length_)
  print(length_)
  x <- paste0('bedtools sort -i ', path.to.gr, ' | bedtools merge -i stdin | bedtools intersect -a stdin -b ', path.to.mut, ' -wa -wb')
  mut.in.df <- tryCatch(read.delim(pipe(x), header = F), error = function(e) NULL)
  if(is.null(mut.in.df)){
    count.mut <- 0
  } else{
    mut.in.df$mut.id <- str_split_fixed(mut.in.df$V7, '_', Inf)[, 1]
    count.mut <- length(unique(mut.in.df$mut.id))  
  }
  print(count.mut)
  return(c(length_, count.mut))
}

TestMutEnrichment <- function(path.to.sample, path.to.population, path.to.mut, subtract = 'False'){
  x <- paste0('awk -F"\\t" \'{print $9,$10,$11,$1"_"$2"_"$15"_"$16"_"$17}\' OFS="\\t" ', path.to.mut, ' | sed \'s/^/chr/\' | tail -n +2 |sort -u | bedtools sort -i stdin > tmp_mut.bed')
  system(x)
  add <- ""
  out <- CountMutations(path.to.sample, 'tmp_mut.bed', add)
  length1 <- out[1]
  count.mut1 <- out[2]
  
  if(subtract == 'True'){
    add <- paste0("| bedtools subtract -a stdin -b ", path.to.sample)
  }
  out <- CountMutations(path.to.population, 'tmp_mut.bed', add)
  length2 <- out[1]
  count.mut2 <- out[2]
  
  system('rm tmp_mut.bed')
  
  p.value <- poisson.test(x = c(count.mut1, count.mut2), T = c(length1, length2), alternative = 't', r = 1)$p.value
  return(list('proc1' = count.mut1 / length1,
              'proc2' = count.mut2 / length2,
              'p.value' = p.value))
}

print(TestMutEnrichment(path.to.sample, path.to.population, path.to.mut, subtract))

