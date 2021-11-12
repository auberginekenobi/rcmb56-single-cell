#!/bin/bash

# lazy directory navigation
cwd=$(pwd)
src=/home/ochapman/Documents/Mesirov/medullo_ecDNA/amplified-intervals-bed


# Convert the cohort to hg38
cd $src/cohort
for f in $(pwd)/*.hg19.bed ; do 
	outfile=$(basename $f .hg19.bed).hg38.bed
	CrossMap.py bed ../../../../anno/hg19ToHg38.over.chain.gz $f $outfile
done 

rm *.unmap

# aggregate
cat *.hg38.bed | \
 bedtools sort -i - | \
 bedtools merge -d 100000 -i - \
 > $cwd/all-mb-patient-ecDNA.hg38.bed 

# Convert the pdx to hg38
cd $src/pdx
for f in $(pwd)/*.hg19.bed ; do
	outfile=$(basename $f .hg19.bed).hg38.bed
	CrossMap.py bed ../../../../anno/hg19ToHg38.over.chain.gz $f $outfile
done

rm *.unmap

# aggregate
cat *.hg38.bed | \
 bedtools sort -i - | \
 bedtools merge -d 100000 -i - \
 > $cwd/all-mb-model-ecDNA.hg38.bed

cd $cwd
cat all-mb-patient-ecDNA.hg38.bed all-mb-model-ecDNA.hg38.bed > all-mb-ecDNA.hg38.bed
