#!/bin/bash
# Make mutations bed file
awk -F"\t" '{print $9,$10,$11,$1"_"$2"_"$15"_"$16"_"$17}' OFS="\t" $2 | sed 's/^/chr/' | tail -n +2 |sort -u | bedtools sort -i stdin > ../out/mutations_$3.bed
# Intersect with left handles
bedtools sort -i $1 | bedtools intersect -a stdin -b ../out/mutations_$3.bed -wa -wb > ../out/panhandles_with_mutations_left_handles.bed
# swap handles
awk -F"\t" '{print $1,$12,$13,$4,$5,$6,$7,$8,$9,$10,$11,$2,$3}' OFS="\t" $1 > ../out/tmp
# Interscet with right handles
bedtools sort -i ../out/tmp | bedtools intersect -a stdin -b ../out/mutations_$3.bed -wa -wb > ../out/panhandles_with_mutations_right_handles.bed
rm ../out/tmp