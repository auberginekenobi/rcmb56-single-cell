#!/bin/bash

cat /home/ochapman/Documents/Mesirov/medullo_ecDNA/amplified-intervals-bed/cohort/*.bed \
/home/ochapman/Documents/Mesirov/medullo_ecDNA/amplified-intervals-bed/pdx/*.bed | \
bedtools sort -i - | \
bedtools merge -i - > \
all-mb-ecDNA-regions.hg38.bed

BEDLENGTH=/home/ochapman/Documents/Mesirov/Software/oscutils/bed-length.py
python $BEDLENGTH all-mb-ecDNA-regions.hg38.bed
