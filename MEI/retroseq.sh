#!/bin/bash

#SBATCH -J retroseq
#SBATCH -o /hpcfs/users/%u/log/retroseq-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants" #"$(dirname "$(readlink -f "$0")")"
modList=("Singularity/3.10.5" "BCFtools/1.17-GCC-11.2.0" "HTSlib/1.17-GCC-11.2.0")

usage()
{
echo "# Script for identifying retrotransposed elements in an aligned Genome using RetroSeq pipeline
# A Genome usually requires 16 cores and time 1 day, mem 8GB
#
# Usage sbatch --array 0-(# of bam files minus 1) $0 -b /path/to/list.of.bamfiles [-a percent-align-identity -o /path/to/output -g /path/to/Genome.fa -r /path/to/repeatsDir ] | [ - h | --help ]
#
# Options
# -b    REQUIRED. List of bam files to call, in a text file with the full path to the bam file provided.
# -a	OPTIONAL. Equivalent to RetroSeq -id option (NOTE: program default is 90 but default for this script is 80)
# -o    OPTIONAL. Path to where you want to find your file output (if not specified /cratchdata1/users/${USER}/retroseq is used)
# -c    OPTIONAL. Path to a RetroSeq config file (if not specified hs38DH.retroseq.phoenix.cfg is used)
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
#
# Original: Atma Ivancevic, 20/08/2018
# Modified: (Date; Name; Description)
# 27/04/2019; Mark Corbett <mark dot corbett at adelaide.edu.au>; Added usage and script flags; centralised perl script and default locations of repeats directory
# 30/05/2025; Mark Corbett <mark dot corbett at adelaide.edu.au>; Updated to use singularity
#
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -b )                    shift
                                        bamList=$1
                                        ;;
                -a )                    shift
                                        pctAlignId=$1
                                        ;;
                -o )                    shift
                                        outputDir=$1
                                        ;;
                -c )                    shift
                                        Config=$1
                                        ;;
                -h | --help )           usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

if [ -z "${bamList}" ]; then # If no list of BAM / CRAM files was provided then do not proceed
    usage
    echo "## ERROR: You need to provide a text file with list of BAM or CRAM files with the full path to the file"
    exit 1
fi

if [ -z "${Config}" ]; then # If no config file specified use the default
    Config=${scriptDir}/configs/hs38DH.retroseq.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi

source ${Config}

bamFile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${bamList}) # Get the BAM file from the list
inputDir=$(dirname "${bamFile}") # Get the input directory from the BAM file path
outPrefix=$(basename "${bamFile}" | sed 's/\.[^.]*$//') # Get the output prefix from the BAM file name. NOTE: This pattern may not match all BAM file names.

if [ -z "${pctAlignId}" ]; then # If no alignment minimum match percentage specified use the default
    pctAlignId=80
fi
if [ -z "${outputDir}" ]; then # If no output directory then use the default
    outputDir=/hpcfs/users/${USER}/retroseq/${outPrefix}
    echo "## INFO: No output directory was specified so the default ${outputDir} will be used."
fi

if [ ! -d "${outputDir}" ]; then
        mkdir -p ${outputDir}
fi
if [ ! -d "${outputDir}/TEdiscovery" ]; then
        mkdir -p ${outputDir}/TEdiscovery
        mkdir -p ${outputDir}/TEcalling
fi

# Do a little summary for the log file
echo "
# Running Retroseq with the following parameters
# File: ${bamFile} 
# Outputs can be found here: ${outputDir}
# Minimum match percentage: ${pctAlignId} %
# Genome: ${Build} (If you have problems with the output make sure this reference matches the BAM file)
# Repeat library: ${repeatsDir}
"
# load modules
for mod in "${modList[@]}"; do
    module load ${mod}
done

# run pipeline

### discovery phase ###
echo "discovering..."

singularity exec --bind ${neuroDir},${outputDir} ${progDir}/${progName} -discover \
-bam ${bamFile} \
-output ${outputDir}/TEdiscovery/${outPrefix}.candidates.tab \
-eref ${repeatsDir}/eref_types.tab \
-align -id ${pctAlignId} > ${outputDir}/${outPrefix}.log 2>&1
echo "discovery complete"

### filtering ###
# filter out candidates near contig starts (likely contamination)
echo "filtering contigs..."
awk -F"\t" '{if ($2>1000) print}' ${outputDir}/TEdiscovery/${outPrefix}.candidates.tab \
> ${outputDir}/TEdiscovery/${outPrefix}.candidates.tab.filtered
echo "filtering complete"

### calling phase ###
echo "calling..."

singularity exec --bind ${neuroDir},${outputDir} ${progDir}/${progName} -call \
-bam ${bamFile} \
-input ${outputDir}/TEdiscovery/${outPrefix}.candidates.tab.filtered \
-ref ${refPath}/${Genome} \
-output ${outputDir}/TEcalling/${outPrefix}.vcf >> ${outputDir}/${outPrefix}.log 2>&1
echo "calling complete"

# compress and index the VCF file
bgzip ${outputDir}/TEcalling/${outPrefix}.vcf && tabix ${outputDir}/TEcalling/${outPrefix}.vcf.gz
