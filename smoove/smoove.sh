#!/bin/bash

## Script for running Smoove which incorporates Lumpy for SVs
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genome/Bams4smoove sbatch smoove.sh

#SBATCH -A robinson
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load Singularity/3.7.1

INDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/alignments/Illumina/genome/Bams4Smoove/Trio
SMOOVEEXE=/hpcfs/groups/phoenix-hpc-neurogenetics/executables
GENOMEDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq
bed=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GFFandBEDfiles/exclude.cnvnator_100bp.GRCh38.20170403.bed
OUTDIR=/hpcfs/users/$USER/outputs/results-smoove/Trio/

# run the thing
echo $(date +"[%b %d %H:%M:%S] Go to dir")
export TMPDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/$USER/tmp/
cd $INDIR

singularity exec $SMOOVEEXE/smoove_v0.2.7.sif smoove call -x --name TRIO --outdir $OUTDIR --exclude $bed --fasta $GENOMEDIR/hs38DH.fa -p 8 --genotype $INDIR/*.bam

echo $(date +"[%b %d %H:%M:%S] All done!")
