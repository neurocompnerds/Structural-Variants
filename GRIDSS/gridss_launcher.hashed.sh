#!/bin/bash
# This is the master script that coordinates job submission for the neurogenetics GRIDSS structural variant pipeline.
## Set hard-coded paths and define functions ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
logDir="/hpcfs/users/${USER}/log"

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

usage()
{
echo "# This is the master script that coordinates job submission for analysis of one or many Illumina whole genomes using GRIDSS.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: An aligned BAM file, samtools, Java 
#
# Usage $0 -p file_prefix -i /path/to/input/bam-file [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX*.bam
# -i	REQUIRED. Path to the BAM file
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/users/${USER}/GRIDSS/\${outPrefix} is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 06/07/2021 
# Modified: (Date; Name; Description)
# 
# 
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -c )            shift
                        Config=$1
                        ;;
        -p )            shift
                        outPrefix=$1
                        ;;
        -i )            shift
                        inputDir=$1
                        ;;
        -o )            shift
                        workDir=$1
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
if [ -z "$Config" ]; then # If no config file specified use the default
    Config=$scriptDir/configs/hs38DH.SV_GRIDSS.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "${outPrefix}" ]; then # If no file prefix specified then do not proceed
    usage
    echo "## ERROR: You need to specify a file prefix (PREFIX) referring to your BAM files eg. PREFIX*.bam."
    exit 1
fi
if [ -z "$inputDir" ]; then # If path to bam file not specified then do not proceed
	usage
	echo "## ERROR: You need to specify the path to your bam files"
	exit 1
fi
# Locate the bam
bamFile=$(find $inputDir/*.bam | grep $outPrefix)
if [ ! -f "$bamFile" ]; then
    echo "## ERROR: BAM file not found in $inputDir"
    exit 1
fi 
if [ -z "$workDir" ]; then # If no output directory then set and create a default directory
	workDir=/hpcfs/users/${USER}/GRIDSS/$outPrefix
	echo "## INFO: Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi

## Launch the job chain ##
#preprocess_job=`sbatch --export=ALL $scriptDir/GRIDSS/gridss_preprocess.sh -c $Config -p $outPrefix -i $inputDir -o $workDir`
#preprocess_job=$(echo $preprocess_job | cut -d" " -f4)
assembly_job=`sbatch --array=0-31 --export=ALL $scriptDir/GRIDSS/gridss_assembly.sh -c $Config -p $outPrefix -i $inputDir -o $workDir`
assembly_job=$(echo $assembly_job | cut -d" " -f4)
sbatch --export=ALL --dependency=afterok:${assembly_job} $scriptDir/GRIDSS/gridss_call.sh -c $Config -p $outPrefix -i $inputDir -o $workDir
