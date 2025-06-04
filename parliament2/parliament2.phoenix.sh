#!/bin/bash
#SBATCH -J CNVByOwl
#SBATCH -o /hpcfs/users/%u/log/parliament2-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00:00
#SBATCH --mem=128GB

# Notification Configuration 
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# A script to call CNV from short read genome squencing

## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants" #"$(dirname "$(readlink -f "$0")")"
modList=("Singularity/3.10.5" "SAMtools/1.17-GCC-11.2.0")

usage()
{
echo "# This script takes BAM files as input and calls structural variants with the Parliament2 package.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: An aligned BAM (or CRAM) file (or files), Singularity and the Parliament2 software.
# NOTE: CRAM or BAM files must be placed in a subdirectory of /hpcfs/groups/phoenix-hpc-neurogenetics
# NOTE: Parliament2 is not respectful of your source data (it will delete your BAM/CRAMs) so the script will copy your files first then try to clean up after itself.
# NOTE: This script is only designed to work with hs38DH genome build so far.  If you have CRAMs and they are not mapped to that build then this won't work.
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

## Load modules ##
for mod in "${modList[@]}"; do
    module load ${mod}
done

readarray -t bamFile < ${bamList} # Read the BAM file list into an array
inputDir=$(dirname "${bamFile[SLURM_ARRAY_TASK_ID]}") # Get the input directory from the BAM file path
outPrefix=$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}" | sed 's/\.[^.]*$//') # Get the output prefix from the BAM file name. NOTE: This pattern may not match all BAM file names.
extn=${bamFile[SLURM_ARRAY_TASK_ID]##*.}

if [ -z "${outputDir}" ]; then # If no output directory then set a default directory
	outputDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/SV/Parliament2/${Build}/${outPrefix}
	echo "## INFO: Using ${outputDir} as the output directory"
else
    outputDir=${outputDir}/${Build}/${outPrefix}
    echo "## INFO: Using ${outputDir} as the output directory for ${bamFile[SLURM_ARRAY_TASK_ID]}"
fi

# Ensure required directories exist
mkdir -p ${outputDir}/in
mkdir -p ${outputDir}/out
cp ${neuroDir}/RefSeq/${Genome} ${outputDir}/in/
cp ${neuroDir}/RefSeq/${Genome}.fai ${outputDir}/in/
refPath=/home/dnanexus/in # Set the reference path to the container's input directory

# Parliament claims to handle CRAMs but in reality it doesn't.
case ${extn} in
    "bam"  )
        baiFile=$(find ${inputDir}/*.bai | grep -w ${outPrefix})
        cp "${bamFile[SLURM_ARRAY_TASK_ID]}" "${outputDir}/in/$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}")"
        cp "${baiFile}" "${outputDir}/in/$(basename "${baiFile}")"
        BF="/home/dnanexus/in/$(basename "${bamFile[SLURM_ARRAY_TASK_ID]}")"
        IndexFile="/home/dnanexus/in/$(basename "${baiFile}")"
        ;;
    "cram" )
        samtools view -T ${neuroDir}/RefSeq/${Genome} -b -@8 -o ${outputDir}/in/${outPrefix}.bam ${bamFile[SLURM_ARRAY_TASK_ID]}
        samtools index ${outputDir}/in/${outPrefix}.bam
        baiFile=$(find ${outputDir}/in/*.bai | grep -w ${outPrefix})
        BF="/home/dnanexus/in/${outPrefix}.bam"
        IndexFile="/home/dnanexus/in/$(basename "${baiFile}")"
        ;;
    * )
        echo "## ERROR: The file ${bamFile[SLURM_ARRAY_TASK_ID]} is not a BAM or CRAM file."
        exit 1
        ;;
esac

#Switch to the output directory because this software dumps all manner of crap in the current working directory if you don't
cd ${outputDir}

singularity exec --bind ${outputDir}/in:/home/dnanexus/in,${outputDir}/out:/home/dnanexus/out \
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
rm -r ${outputDir}/in
rm ${outputDir}/ref.fa ${outputDir}/ref.fa.fai # because what I really want is multiple copies of the reference genome in every directory!
