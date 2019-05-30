#!/bin/bash
awk -F"\t" '{print $1,$11,$12,$4,$5,$6,$7,$8,$9,$10,$2,$3}' OFS="\t" $1 > alignments_whole_human_genome_processed_bed12_right.bed
