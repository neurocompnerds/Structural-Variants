#!/bin/bash

#SBATCH -J mergeDedupBAM
#SBATCH -o /hpcfs/users/%u/log/mergeDedupBAM-slurm-%j.out
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
modList=("HTSlib/1.17-GCC-11.2.0" "SAMtools/1.17-GCC-11.2.0")

usage()
{
echo "# Script for sorting and marking duplicate reads in BAM files
# Requires: samtools, sambamba
#
# Usage sbatch --array 0-\"number of samples -1\" $0 -i sample.file.txt [ -c /path/to/Config.cfg -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields: SampleName.
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
                    sampleFile=$1
                    ;;
        -c )        shift
                    enviroCfg=$1
                    ;;
        -o )        shift
                    outDir=$1
                    ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )         usage
                    exit 1
    esac
    shift
done

## Pre-flight checks ##
if [ -z "${sampleFile}" ]; then # If no sample file specified, display usage and exit
    echo "## ERROR: No sample file specified."
    usage
    exit 1
fi

readarray -t Sample <<< $(grep -v "^#" ${sampleFile} | cut -f1) 

if [ -z "${enviroCfg}" ]; then # If there is no config file then use the default.
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,wisecondorx,configs,g')"
    enviroCfg="${configDir}/hs38DH.SV_WISECONDORX.phoenix.cfg"
    echo "## INFO: the default config file ${enviroCfg} has been selected."
fi

source ${enviroCfg}

tmpDir=${baseTmpDir}/${Sample[SLURM_ARRAY_TASK_ID]}
if [ ! -d "${tmpDir}" ]; then
    echo "## ERROR: The temporary directory ${tmpDir} does not exist.  This should have been created in the mapping step.  Please check your log files for errors."
    exit 1
fi

## Create essential directories ##
if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
    echo "## INFO: ${outDir} has been created as the output directory"
fi
for dir in qc align; do
    if [ ! -d "${outDir}/${dir}" ]; then
        mkdir -p ${outDir}/${dir}
        echo "## INFO: Created ${outDir}/${dir}"
    fi
done

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Mark duplicates ##
${sambambaProg} markdup \
    -t 24 \
    -l 5 \
    --tmpdir=${tmpDir} \
    --overflow-list-size 1000000 \
    --hash-table-size 1000000 \
    ${tmpDir}/${Sample[SLURM_ARRAY_TASK_ID]}.*.samsort.bwa.${Build}.bam \
    ${outDir}/align/${Sample[SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.${Build}.bam

## Clean up space ##
if [ -f "${outDir}/align/${Sample[SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.${Build}.bam" ]; then
    rm ${tmpDir}/${Sample[SLURM_ARRAY_TASK_ID]}.*.samsort.bwa.${Build}.bam
else
    echo "## ERROR: Duplicate marking or earlier stage failed!"
    exit 1
fi

## Collect stats ##
echo "# Flagstats" > ${outDir}/qc/${Sample[SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
samtools flagstat ${outDir}/align/${Sample[SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.${Build}.bam >> ${outDir}/qc/${Sample[SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
