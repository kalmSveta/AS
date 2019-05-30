#!/bin/bash

bedtools flank -i $2 -g $3 -b $5 > flanking_regions_$1.bed 
echo 'made flanking regions'
tail -n +6 $4 | awk -F"\t" '($3 == "exon")'| grep 'gene_type \"protein_coding\";' | bedtools subtract -a flanking_regions_$1.bed -b stdin| bedtools sort -i stdin | bedtools map -a stdin -b track.bg -c 4 -o sum -null 0 > flanking_sum_$1.bed 
echo 'made flanking sum'
bedtools sort -i $2 |  bedtools map -a stdin -b track.bg -c 4 -o sum -null 0 > panhandles_sum_$1.bed  
echo 'made panhandles sum'
