#!/bin/bash
# Make DECIPHER common del and dup bed files
wget -N https://decipher.sanger.ac.uk/files/downloads/population_cnv.txt.gz
gunzip population_cnv.txt.gz
awk -F"\t" '$6 > 0.01 && $15 > 200 {print}' population_cnv.txt | grep -v \# | cut -f2,3,4 | sort -k1,1 -k2,2n | bedtools merge -i stdin > common.DECIPHER.dels.gt200.fq0.01.b37.bed
awk -F"\t" '$9 > 0.01 && $15 > 200 {print}' population_cnv.txt | grep -v \# | cut -f2,3,4 | sort -k1,1 -k2,2n | bedtools merge -i stdin > common.DECIPHER.dups.gt200.fq0.01.b37.bed
sed 's/^/chr/g' common.DECIPHER.dels.gt200.fq0.01.b37.bed > common.DECIPHER.dels.gt200.fq0.01.hg19.bed
sed 's/^/chr/g' common.DECIPHER.dups.gt200.fq0.01.b37.bed > common.DECIPHER.dups.gt200.fq0.01.hg19.bed
gzip population_cnv.txt

# dgv for hg19
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz
gunzip dgvMerged.txt.gz
awk -F"\t" '$20/$18 > 0.01 && $18 > 200 {print}' dgvMerged.txt | grep -v \# | cut -f2,3,4 | sort -k1,1 -k2,2n | bedtools merge -i stdin > common.DGV.dels.gt200.fq0.01.hg19.bed
awk -F"\t" '$19/$18 > 0.01 && $18 > 200 {print}' dgvMerged.txt | grep -v \# | cut -f2,3,4 | sort -k1,1 -k2,2n | bedtools merge -i stdin > common.DGV.dups.gt200.fq0.01.hg19.bed
sed 's/^chr//g' common.DGV.dels.gt200.fq0.01.hg19.bed > common.DGV.dels.gt200.fq0.01.b37.bed
sed 's/^chr//g' common.DGV.dups.gt200.fq0.01.hg19.bed > common.DGV.dups.gt200.fq0.01.b37.bed
gzip dgvMerged.txt
