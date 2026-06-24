#!/bin/bash -l
#SBATCH -J bam2npz
#SBATCH -o /hpcfs/users/%u/log/WC-convert-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=06:00:00
#SBATCH --mem=128GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for converting bam files to npz arrays for analysis with WisecondorX
## Set hard coded paths and variables ##

modList=("Anaconda3/2024.06-1")
usage()
{
echo "# Script for converting bam files to npz arrays for analysis with WisecondorX
# Requires: WisecondorX conda environment please see the readme for instructions on how to set this up.
#
# Usage sbatch --array 0-\"number of samples -1\" -i sample.file.txt [ -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields: SampleName.
# -c    OPTIONAL. Path to the configuration file.
# -o    OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
# -h or --help    Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett; 18/06/2026
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

readarray -t Sample <<< $(grep -v "^#" ${sampleFile} | cut -f1)

if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi
if [ ! -d "${outDir}/align" ]; then
    echo "## ERROR: ${outDir}/align not found. The script expects to find your BAM files at this location."
    exit 1
fi
if [ ! -d "${outDir}/wc" ]; then
    mkdir -p ${outDir}/wc
    echo "## INFO: Created WisecondorX output directory: ${outDir}/wc"
fi

## Create npz arrays from the BAM files ##
conda activate WisecondorX
wisecondorx convert \
    ${outDir}/align/${Sample[SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.${Build}.bam \
    ${outDir}/wc/${Sample[SLURM_ARRAY_TASK_ID]}.merge.dedup.npz \
    --threads ${SLURM_NTASKS}
conda deactivate
