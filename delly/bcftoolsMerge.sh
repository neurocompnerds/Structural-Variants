#!/bin/bash

## Script for merging genotyped SVs with bcftools
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/hpcfs/users/$USER/outputs/SVcalling sbatch bcftoolsMerge.sh

#SBATCH -A robinson
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load BCFtools/1.9-foss-2016b

INDIR=/hpcfs/users/$USER/outputs/SVcalling

# run the thing
echo $(date +"[%b %d %H:%M:%S] Go to dir")
cd $INDIR
pwd

echo $(date +"[%b %d %H:%M:%S] Set bcf list")
genotypedList=$(find *.geno.bcf)
echo $genotypedList

echo $(date +"[%b %d %H:%M:%S] Merge all genotyped bcfs")
bcftools merge -m id -O b -o merged.bcf $genotypedList

echo $(date +"[%b %d %H:%M:%S] All done!")
