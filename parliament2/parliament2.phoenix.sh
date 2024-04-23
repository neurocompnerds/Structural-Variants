#!/bin/bash
#SBATCH -J OwlsCallCNV
#SBATCH -o /hpcfs/users/%u/log/parliament2-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=1-00:00:00
#SBATCH --mem=16GB

# Notification Configuration 
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# A script to call CNV from short read genome squencing
## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
variantDir="/hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2"
modList=("arch/skylake" "Singularity/3.7.4")

usage()
{
echo "# This script takes BAM files as input and calls structural variants with the Parliament2 package.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: An aligned BAM file (or files), Singularity
#
# Usage sbatch --array 0-(n-1 bam files) $0 -i /path/to/input/folder_of_bams [ -b listOfbamFiles -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -i	REQUIRED. Path to the folder containing the BAM files
# -b    OPTIONAL. A list of bam files to call, if not included all of the bam files in the folder will be called
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2/\${BUILD}/\${outPrefix} is used)
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
        -b )            shift
                        outPrefix=$1
                        ;;
        -i )            shift
                        inputDir=$1
                        ;;
        -o )            shift
                        outputDir=$1
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
bamFile=$(find $inputDir/*.bam | grep -w $outPrefix)
if [ ! -f "$bamFile" ]; then
    echo "## ERROR: BAM file not found in $inputDir"
    exit 1
fi 
baiFile=$(find $inputDir/*.bai | grep -w $outPrefix)
if [ ! -f "$baiFile" ]; then
    echo "## ERROR: The BAM index file was not found in $inputDir"
    exit 1
fi 
if [ -z "$outputDir" ]; then # If no output directory then set and create a default directory
	outputDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2/${BUILD}/
	echo "## INFO: Using $outputDir as the output directory"
fi
if [ ! -d "$outputDir" ]; then
    mkdir -p $outputDir
fi
tmpDir=$baseTmpDir/$outPrefix
if [ ! -d "$tmpDir" ]; then
    mkdir -p $tmpDir
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

singularity run -v ${inputDir}:/home/dnanexus/in -v ${outputDir}:/home/dnanexus/out dnanexus/parliament2:<TAG> --bam ${bamFile} --bai ${baiFile} --fai <REFERENCE_INDEX> -r <REFERENCE_NAME> <OPTIONAL_ARGUMENTS>