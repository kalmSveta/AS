#!/bin/bash
length=$3
# select intronic regions
if [[ "$4" == "True" ]]; then
    python ../lib/gtftools.py -d ../data/introns_gtftools.bed $1 -c $5
    awk -F"\t" '{print $1,$2,$3,$4}' OFS="\t" ../data/introns_gtftools.bed > ../data/tmp
    mv ../data/tmp ../data/introns_gtftools.bed
    echo "selected always intronic regions"
else
    python ../lib/gtftools.py -i ../data/introns_gtftools.bed $1 -c $5
    awk -F"\t" '{print $1,$2,$3,$5}' OFS="\t" ../data/introns_gtftools.bed > ../data/tmp
    mv ../data/tmp ../data/introns_gtftools.bed
    echo "selected introns"
fi
# add chr, make0-based
sed 's/^/chr/' ../data/introns_gtftools.bed | awk -F"\t" '{$2 = $2 - 1; print}' OFS="\t" | awk -F"\t" '{$3 = $3 - 1; print}' OFS="\t" > ../data/tmp
mv ../data/tmp ../data/introns_gtftools.bed

# intersect with conserved regions, merge, intersect with coding genes
awk -F"\t" '{print $2,$3,$4}' OFS="\t" $2 | bedtools intersect -a ../data/introns_gtftools.bed -b stdin |  awk -F"\t" '{print $1"_"$4,$2,$3,$4}' OFS="\t" | bedtools sort -i stdin|bedtools merge -c 4 -o distinct -i stdin| awk -F"\t" '{n = split($1, a, "_"); print a[1],$2,$3,$4}' OFS="\t" | bedtools intersect -a stdin -b ../data/coding_genes.bed -wa -wb | awk -F"\t" '{if ($4 == $8) {print $1,$2,$3,$11,$6,$7}}' OFS="\t"  > ../data/tmp

# select only long enough
awk -v len="$length" -F"\t" '{if($3-$2 >= len) { print }}' OFS="\t" ../data/tmp > ../data/conin_gtftools_long_coding.tsv
# add header
echo -e "chr\tstart_interval\tend_interval\tstrand\tstart_gene\tend_gene" | cat - ../data/conin_gtftools_long_coding.tsv > ../data/conin_gtftools_long_coding_final.tsv
rm ../data/tmp
rm ../data/conin_gtftools_long_coding.tsv
