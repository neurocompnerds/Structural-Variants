#!/bin/bash

#SBATCH -J goamogogoamo
#SBATCH -o /hpcfs/users/%u/log/mosdepth-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=06:00:00
#SBATCH --mem=32GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for calculating read depth from genome sequencing data
## Set hard coded paths and variables ##

usage()
{
echo "# Script for calculating read depth from genome sequencing data
# Requires: mosdepth
#
# Usage sbatch --array 0-\"number of samples -1\" $0 -i sample.file.txt [ -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required field: SampleName.
# -c    OPTIONAL. Path to the configuration file.
# -o    OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
# -h or --help    Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: By Mark Corbett, 18/06/2026
# Modified: (Date; Name; Description)
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
if [ -z "${enviroCfg}" ]; then # If there is no config file then use the default.
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,wisecondorx,configs,g')"
    enviroCfg="${configDir}/hs38DH.SV_WISECONDORX.phoenix.cfg"
    echo "## INFO: the default config file ${enviroCfg} has been selected."
fi

source ${enviroCfg}

if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi

## Parse the input file ##
Sample=$((grep -v "^#" ${sampleFile} | cut -f1)) 

## Create essential directories ##

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
    echo "## INFO: ${outDir} has been created as the output directory"
fi
if [ ! -d "${outDir}/qc" ]; then
    mkdir -p ${outDir}/qc
fi

## run mosdepth ##
${mosdepthProg} ${outDir}/qc/${Sample[SLURM_ARRAY_TASK_ID]}.depth -n  -x -t 4 ${outDir}/${Sample[SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.${Build}.bam
