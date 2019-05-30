#!/bin/bash
path=${1//$2/}
echo ../out/handles_and_mut/"$path"_out.txt
vcftools --gzvcf $1 --recode --stdout --remove-indels | bedtools intersect -a $3 -b stdin -wa -wb  > ../out/handles_and_mut/"$path"_out.txt
for file in ../out/handles_and_mut/*out.txt; do gzip $file; done
