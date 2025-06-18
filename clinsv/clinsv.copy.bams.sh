#!/bin/bash
#SBATCH -J clinSVcp
#SBATCH -o /hpcfs/users/%u/log/clinSVcp-slurm-%j.out

#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=5:30:00
#SBATCH --mem=36GB

# Notification Configuration 
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# This is a helper script that sets up files for the ClinSV structural variant pipeline.
## Set hard-coded paths and define functions ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
modList=("SAMtools/1.17-GCC-11.2.0")

usage()
{
echo "# This is a helper script that coordinates data ingress for analysis of one or many Illumina whole genomes using ClinSV.
# The script should be launched using the clinsv_launcher.sh script, but you can run it independently if needed.  
# Requires: A short list of aligned CRAM or BAM file 
# Note: All output is written to /hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/test_run
#
# Usage sbatch --array 0-(number of bam files - 1) $0 -b /path/to/input/bam-file-list.txt [ -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -b	REQUIRED. Path to the file containing list of CRAM/BAM files (with their full path) to be processed.  The file should contain one CRAM/BAM file per line.
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 06/07/2021 
# Modified: (Date; Name; Description)
# Note: Refer to github for edit history
# 
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -c )            shift
                        Config=$1
                        ;;
        -b )            shift
                        bamList=$1
                        ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )             usage
                        exit 1
    esac
    shift
done

## Pre-flight checks ##
if [ -z "${bamList}" ]; then # If no list of BAM / CRAM files was provided then do not proceed
    usage
    echo "## ERROR: You need to provide a text file with list of BAM or CRAM files with the full path to the file"
    exit 1
fi

if [ -z "${Config}" ]; then # If no config file specified use the default
    Config="${scriptDir}/configs/hs38DH.SV_ClinSV.phoenix.cfg"
    echo "## INFO: Using the default config ${Config}"
fi

source "${Config}"

## Load modules ##
for mod in "${modList[@]}"; do
    module load ${mod}
done

# Convert CRAM to BAM and copy files over to the clivar directory
readarray -t bamFile < "${bamList}"
outPrefix=$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}" | sed 's/\.[^.]*$//') # Remove the file extension to get the output prefix
inputDir=$(dirname "${bamFile[SLURM_ARRAY_TASK_ID]}")
extn="${bamFile[SLURM_ARRAY_TASK_ID]##*.}"
case ${extn} in
    "bam"  )
        if [ -f "${neuroDir}/clinsv/${outPrefix}.bam" ]; then
            echo "## INFO: The file ${neuroDir}/clinsv/${outPrefix}.bam already exists so I'm going to use that to save some time."
            exit 0
        fi
        baiFile=$(find "${inputDir}/*.bai" | grep -w "${outPrefix}")
        cp "${bamFile[SLURM_ARRAY_TASK_ID]}" "${neuroDir}/clinsv/"
        cp "${baiFile}" "${neuroDir}/clinsv/"
        ;;
    "cram" )
        if [ -f "${neuroDir}/clinsv/${outPrefix}.bam" ]; then
            echo "## INFO: The file ${neuroDir}/clinsv/${outPrefix}.bam already exists so I'm going to use that to save some time."
            exit 0
        fi        
        samtools view -T "${neuroDir}/RefSeq/${Genome}" -b -@8 -o "${neuroDir}/clinsv/${outPrefix}.bam" "${bamFile[SLURM_ARRAY_TASK_ID]}"
        samtools index "${neuroDir}/clinsv/${outPrefix}.bam"
        ;;
    * )
        echo "## ERROR: The file ${bamFile[SLURM_ARRAY_TASK_ID]} is not a BAM or CRAM file."
        exit 1
        ;;
esac
