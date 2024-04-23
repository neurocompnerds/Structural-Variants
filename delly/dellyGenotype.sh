#!/bin/bash

## Script for genotyping SVs with delly
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genome/allGenomes SITELIST=/hpcfs/users/$USER/outputs/SVcalling/dellyOut/sites.bcf sbatch --array 0-19 dellyGenotype.sh

#SBATCH -p icelake,a100cpu
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

INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genome/AGRF_KCCG_WGS/short_bams
SITELIST=/hpcfs/users/$USER/outputs/SVcalling/sites.bcf
DELLYEXE=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/delly/delly
OUTDIR=$FASTDIR/outputs
GENOMEDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# run the thing

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting delly genotyping")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

$DELLYEXE call \
-g ${GENOMEDIR}/hs38DH.fa \
-v ${SITELIST} \
-o ${OUTDIR}/SVcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.geno.bcf \
-x ${GENOMEDIR}/human.hg38.excl.tsv \
${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
2> ${FASTDIR}/slurmLOG/dellyCalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo $(date +"[%b %d %H:%M:%S] All done!")
