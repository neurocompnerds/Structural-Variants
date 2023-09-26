#!/bin/bash

## Script for calling SVs with delly
## Date: 28 March 2018 
##
## Example usage:
## INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genomes/ID_Jan2021 sbatch --array 0-3 dellyCall.sh

#SBATCH -A robinson
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3-00:00
#SBATCH --mem=32GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/hpcfs/users/$USER

DELLYEXE=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/delly/delly
OUTDIR=$FASTDIR/outputs
GENOMEDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/
INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genome/ID_Jan2021

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# load modules
module load arch/haswell
module load BCFtools/1.9-foss-2016b

# run the thing

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting delly")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

$DELLYEXE call \
-g ${GENOMEDIR}/hs38DH.fa \
-o ${OUTDIR}/SVcalling/dellyOut/${QUERIES[$SLURM_ARRAY_TASK_ID]}.bcf \
-x ${GENOMEDIR}/human.hg38.excl.tsv \
${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
2> ${FASTDIR}/slurmLOG/dellyCalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo $(date +"[%b %d %H:%M:%S] All done!")
