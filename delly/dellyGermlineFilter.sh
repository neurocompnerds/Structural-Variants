#!/bin/bash

## Script for filtering merged SVs with the delly germline filter
## Requires at least 20 unrelated samples
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/hpcfs/users/$USER/outputs/SVcalling sbatch dellyGermlineFilter.sh

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

# set variables
INDIR=/hpcfs/users/$USER/outputs/SVcalling
DELLYEXE=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/delly/delly

# load modules
module load arch/haswell
module load BCFtools/1.9-foss-2016b

# run the thing
echo $(date +"[%b %d %H:%M:%S] Go to dir")
cd $INDIR
pwd

echo $(date +"[%b %d %H:%M:%S] Check for merged bcf")
ls merged.bcf

echo $(date +"[%b %d %H:%M:%S] Index bcf")
bcftools index merged.bcf

echo $(date +"[%b %d %H:%M:%S] Apply delly germline filter")
$DELLYEXE filter -f germline -o germline.bcf merged.bcf
# if you want to keep only variants where FILTER=PASS
# add -p to the above command

echo $(date +"[%b %d %H:%M:%S] Also output vcf")
bcftools view germline.bcf > germline.vcf

echo $(date +"[%b %d %H:%M:%S] All done!")
