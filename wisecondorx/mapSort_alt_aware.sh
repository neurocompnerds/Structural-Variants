#!/bin/bash

#SBATCH -J MapSort
#SBATCH -o /hpcfs/users/%u/log/mapSort-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --time=06:00:00
#SBATCH --mem=148GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for mapping Illumina pair-end sequence data
## Set hard coded paths and variables ##
modList=("BWA/0.7.17-GCCcore-11.2.0" "HTSlib/1.17-GCC-11.2.0" "SAMtools/1.17-GCC-11.2.0")
usage()
{
echo "# Script for mapping Illumina pair-end sequence data
# Requires: BWA 0.7.x, samtools, sambamba
#
# Usage sbatch --array=0-\"number of read_groups -1\" $0 -i input.file.txt [ -c /path/to/Config.cfg -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields SampleName, /path/to/R1.fastq.gz, /path/to/R2.fastq.gz, sample_type (test or control). Plus optional fields for library, platform and platform unit.
# -c    OPTIONAL. Path to the configuration file.
# -o    OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
# -h or --help    Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 03/06/2026; Mark Corbett; Fork from https://github.com/speleonut/map-n-call/blob/master/GATK4/mapSortDedupMarkIndels_alt_aware_phoenix.sh and converted to an array-capable script.
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -i )        shift
                    inputFile=$1
                    ;;
        -c )        shift
                    enviroCfg=$1
                    ;;
        -o )        shift
                    outDir=$1
                    ;;
        -h | --help )       usage
                            exit 0
                            ;;
        * )         usage
                    exit 1
    esac
    shift
done

## Pre-flight checks ##
if [ -z "${inputFile}" ]; then # If no input file specified, display usage and exit
    usage
    exit 1
fi
if [ -z "${enviroCfg}" ]; then # If there is no config file then use the default.
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,wisecondorx,configs,g')"
    enviroCfg="${configDir}/hs38DH.SV_WISECONDORX.phoenix.cfg"
    echo "## INFO: the default config file ${enviroCfg} has been selected."
fi

source ${enviroCfg}
source ${scriptDir}/utilities/shared.functions.sh

if [ -z "${Build}" ]; then # If Build is not in the config then let the user know and exit
    echo "## ERROR: Genome Build not specified in the config file.  Please check your config file and make sure the Build variable is set to one of the following: hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 or CHM13v2."
    exit 1
else
    test_genome_build # This function is defined in the configs/ file and sets the variable genomeType to either "has_alt_contigs" or "no_alt_contigs" depending on the genome Build specified.
fi
if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi

## Parse the input file ##
# This assumes you use the provided template and keep all of the columns in the right order.

readarray -t Sample <<< $(grep -v "^#" ${inputFile} | cut -f1) 
readarray -t seqFile1 <<< $(grep -v "^#" ${inputFile} | cut -f2) 
readarray -t seqFile2 <<< $(grep -v "^#" ${inputFile} | cut -f3) 
readarray -t LB <<< $(grep -v "^#" ${inputFile} | cut -f6) # Library information (optional)
readarray -t PL <<< $(grep -v "^#" ${inputFile} | cut -f7) # Platform information (optional)
readarray -t PU <<< $(grep -v "^#" ${inputFile} | cut -f8) # Platform unit information (optional)
ID=$(zcat ${seqFile1[SLURM_ARRAY_TASK_ID]} | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}') # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq

## Check data and add defaults if needed
if [ -z "${Sample[SLURM_ARRAY_TASK_ID]}" || "${seqFile1[SLURM_ARRAY_TASK_ID]}" || "${seqFile2[SLURM_ARRAY_TASK_ID]}" ]; then # If the Sample or fastq are not specified that's a hard fail
    echo "## ERROR: Missing Sample name or fastq file.  Please check your input file and ensure these are specified
    Sample: ${Sample[SLURM_ARRAY_TASK_ID]}
    R1 fastq file: ${seqFile1[SLURM_ARRAY_TASK_ID]}
    R2 fastq file: ${seqFile2[SLURM_ARRAY_TASK_ID]}"
    exit 1
fi
if [ -z "${LB[SLURM_ARRAY_TASK_ID]}" ]; then
    LB="LC_WGS"
    echo "## INFO: Setting Library name to ${LB}."
fi
if [ -z "${PL[SLURM_ARRAY_TASK_ID]}" ]; then
    PL="ILLUMINA"
    PU=$(zcat ${seqFile1[SLURM_ARRAY_TASK_ID]} | head -n 1 | awk -F : '{OFS="."; print $3, $4}') # FLOWCELL.LANE
    echo "## INFO: Setting platform to ${PL} and the platform unit to ${PU}. This is a guess but it doesn't affect your results."
fi

## Create essential directories ##

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
    echo "## INFO: ${outDir} has been created as the output directory"
fi

tmpDir=${baseTmpDir}/${Sample[SLURM_ARRAY_TASK_ID]}
if [ ! -d "${tmpDir}" ]; then
    mkdir -p ${tmpDir}
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
# Map reads to genome using BWA-MEM
# do NOT use the -M option when aligning to GRCh38 as alignment is alt-aware
# -K flag asks bwa-mem to load a fixed number of bases into RAM so enables reproducibility
 
cd ${tmpDir}
bwa mem -K 100000000 -t 24 -R "@RG\tID:${ID}\tLB:${LB[SLURM_ARRAY_TASK_ID]}\tPL:${PL}\tSM:${Sample[SLURM_ARRAY_TASK_ID]}" ${BWAindex} ${seqFile1[SLURM_ARRAY_TASK_ID]} ${seqFile2[SLURM_ARRAY_TASK_ID]} |\
samtools view -bT ${BWAindex} - |\
samtools sort -l 5 -m 4G -@24 -T ${Sample[SLURM_ARRAY_TASK_ID]} -o ${tmpDir}/${Sample[SLURM_ARRAY_TASK_ID]}.${ID}.samsort.bwa.${Build}.bam -
