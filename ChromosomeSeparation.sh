#!/bin/sh

mkdir -p chromosome

for i in `seq 22` X Y

do

  cat refFlat.hg38.txt | grep -w chr${i} > chromosome/refFlat.hg38.chr${i}.txt

done
