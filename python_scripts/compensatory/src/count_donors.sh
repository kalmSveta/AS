#!/bin/bash
#For each file for each row counts number of donors with at least one 1 and donors with 1|1, subtract these to obtain number of donors with mutation in at least one chr   
for file in $1*out.txt.gz; do 
	gunzip $file 
	file_name=${file//txt.gz/txt}
	echo $file_name 
	cut -f 15- $file_name | awk '{print $0,"\t",gsub(/1/,"")}' | awk '{print $0,"\t",gsub(/1[|]1/,"")}'| awk '{print $(NF-1) - $NF}' | paste $file_name - > "$file_name".counts ; 
	gzip $file_name
done
#takes info only about mutations and number of donors
for file in $1*counts; do 
	echo $file
	awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$9,$10,$(NF)}' < $file > $file.2 ; 
done
rm $1/*counts
