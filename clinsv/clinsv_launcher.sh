#!/bin/bash
## Set hard-coded paths and define functions ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
logDir="/hpcfs/users/${USER}/log"

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find the log information from this slurm job here: ${logDir}"
fi

usage()
{
echo "# This is the master script that coordinates job submission for analysis of one or many Illumina whole genomes using ClinSV.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: A short list of aligned CRAM or BAM file 
# Note: All output is written to /hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/test_run
#
# Usage sbatch $0 -b /path/to/input/bam-file-list.txt [ -c /path/to/config.cfg ] | [ - h | --help ]
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
# Check if there is another run in progress. This is not very sophisticated but will work OK with our level of throughput.
if [ -f "/hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/clinsv.lock" ]; then
    echo "## INFO: Sorry it looks like there is already a run in progress. Ask neurogenetics slack / email to find out who is running it and get them to let you know when it is done."
    echo "## INFO: If you are sure that no one is running a job (i.e. you checked) then delete the lock file /hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/.clinsv.lock and try again."
    exit 0
fi
touch /hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/clinsv.lock # Create a lock file to prevent multiple runs at the same time

if [ -z "${bamList}" ]; then # If no list of BAM / CRAM files was provided then do not proceed
    usage
    echo "## ERROR: You need to provide a text file with list of BAM or CRAM files with the full path to the file"
    exit 1
fi
fileCount=$(wc -l < "${bamList}")
jobCount=$((${fileCount} - 1)) # SLURM_ARRAY_TASK_ID starts at 0 so we need to subtract 1 from the file count

if [ -z "${Config}" ]; then # If no config file specified use the default
    Config=${scriptDir}/configs/hs38DH.SV_ClinSV.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi

source ${Config}

## Request some jerbs ##
bamCopyJob=`sbatch --array 0-${jobCount} --export=Config=${Config},bamList=${bamList} ${scriptDir}/clinsv.copy.bams.sh -b ${bamList} -c ${Config} | awk '{print $4}'`
sbatch --dependency=afterok:${bamCopyJob} --export=Config=${Config} ${scriptDir}/clinsv.run.sh -c ${Config}
