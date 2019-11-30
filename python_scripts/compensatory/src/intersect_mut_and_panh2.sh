#!/bin/bash
for mut_file in $1*.vcf.gz; do ../../../tools/vcft/bin/vcftools --gzvcf $mut_file --recode --stdout --remove-indels | bedtools intersect -a $3 -b stdin -wa -wb  > "$mut_file"_out.txt; done

