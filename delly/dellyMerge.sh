#!/bin/bash

## Script for merging SVs with delly
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/hpcfs/users/$USER/outputs/SVcalling sbatch dellyMerge.sh

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/hpcfs/users/$USER

INDIR=/hpcfs/users/$USER/outputs/SVcalling/
DELLYEXE=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/delly/delly
OUTDIR=$FASTDIR/outputs
GENOMEDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/

# load modules
module load arch/haswell
module load BCFtools/1.9-foss-2016b 

# run the thing
echo $(date +"[%b %d %H:%M:%S] Got to dir")
cd $INDIR
pwd

# filter bams to only retain chr1-22,X,Y
# because different hg19 references have diff alt contigs
# and this will make delly merge fail
echo $(date +"[%b %d %H:%M:%S] Restricting variants to chr1-22,X,Y")
for i in *.bam.bcf; do bcftools view $i --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY -O b -o ${i%.bam.bcf}.filtered.bcf; done

# set the list of bams to merge
echo $(date +"[%b %d %H:%M:%S] Set bcf list")
bcfList=$(find *filtered.bcf)
echo $bcfList

# merge bams
# n increases the max SV size to 100000000
echo $(date +"[%b %d %H:%M:%S] Merge all bcfs")
$DELLYEXE merge -n 100000000 -o sites.bcf $bcfList

echo $(date +"[%b %d %H:%M:%S] All done!")
