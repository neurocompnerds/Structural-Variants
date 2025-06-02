#!/bin/bash
#SBATCH -J CNVByOwl
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
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants" #"$(dirname "$(readlink -f "$0")")"
modList=("Singularity/3.10.5")

usage()
{
echo "# This script takes BAM files as input and calls structural variants with the Parliament2 package.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: An aligned BAM (or CRAM) file (or files), Singularity
# NOTE: CRAM or BAM files must be placed in a subdirectory of /hpcfs/groups/phoenix-hpc-neurogenetics
# NOTE: Parliament2 is not respectful of your source data (it will delete your BAM/CRAMs) so the script will copy your files first then try to clean up after itself.
#
# Usage sbatch --array 0-(n-1 bam files) $0 -b listOfbamFiles [-o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -b    REQUIRED. A list of bam files to call, in a text file with the full path to the bam file provided.
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2/\${BUILD}/\${outPrefix} is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 30/05/2025 
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
if [ -z "${bamList}" ]; then # If no list of BAM / CRAM files was provided then do not proceed
    usage
    echo "## ERROR: You need to provide a text file with list of BAM or CRAM files with the full path to the file"
    exit 1
fi

if [ -z "${Config}" ]; then # If no config file specified use the default
    Config=${scriptDir}/configs/hs38DH.Parliament.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi

source ${Config}

readarray -t bamFile < ${bamList} # Read the BAM file list into an array
inputDir=$(dirname "${bamFile[SLURM_ARRAY_TASK_ID]}") # Get the input directory from the BAM file path
outPrefix=$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}" | sed 's/\.[^.]*$//') # Get the output prefix from the BAM file name. NOTE: This pattern may not match all BAM file names.
baiFile=$(find ${inputDir}/*.bai | grep -w ${outPrefix})
if [ ! -f "$baiFile" ]; then
    baiFile=$(find ${inputDir}/*.crai | grep -w ${outPrefix})
    if [ ! -f "${baiFile}" ]; then
        echo "## ERROR: The BAM or CRAM index for ${bamFile[SLURM_ARRAY_TASK_ID]} was not found."
        exit 1
    fi
fi

# Protect your source data
cp "${bamFile[SLURM_ARRAY_TASK_ID]}" "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}")"
cp "${baiFile}" "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${baiFile}")"

# swap in the annoying hard coded dnanexus file path (yes these lines make me feel dirty too).
BF=$(echo "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}")" | sed 's,\/hpcfs\/groups\/phoenix-hpc-neurogenetics,\/home\/dnanexus\/in,g') 
IndexFile=$(echo "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${baiFile}")" | sed 's,\/hpcfs\/groups\/phoenix-hpc-neurogenetics,\/home\/dnanexus\/in,g') 

if [ -z "${outputDir}" ]; then # If no output directory then set a default directory
	outputDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2/${Build}/${outPrefix}
	echo "## INFO: Using ${outputDir} as the output directory"
fi
# Ensure required directories exist
if [ ! -d "${outputDir}" ]; then
    mkdir -p ${outputDir}
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load ${mod}
done

singularity exec --bind ${neuroDir}:/home/dnanexus/in,${outputDir}:/home/dnanexus/out \
    ${progDir}/${progName} \
    --bam ${BF} \
    --bai ${IndexFile} \
    --prefix ${outPrefix} \
    --fai ${refPath}/${Genome}.fai \
    -r ${refPath}/${Genome} \
    --manta --cnvnator --lumpy \
    --delly_deletion --delly_duplication --delly_inversion --delly_insertion \
    --genotype

# Clean up a bit
rm "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}")"
rm "${neuroDir}/alignments/Illumina/genome/bams4parliament2/$(basename "${baiFile}")"
