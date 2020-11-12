#!/bin/sh
set -e

file=$1

prefix=$(basename $file | sed -e 's/.[^.]*$//')
mkdir -p out

#file=refFlat.GRCh38.gene.name.txt
#out_prefix=refFlat.GRCh38.gene.name

cut -f3 $file | sort | uniq -c > out/chromosome_count.txt

for chr in $(cat out/chromosome_count.txt | awk '{print $2}')
do
  cat ${file} | awk '{if($3==$chr){print}}' > out/${prefix}_${chr}.txt
done
