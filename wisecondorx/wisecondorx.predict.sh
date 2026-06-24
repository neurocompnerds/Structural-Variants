#!/bin/bash -l
#SBATCH -J wcPredict
#SBATCH -o /hpcfs/users/%u/log/WC-predict-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=06:00:00
#SBATCH --mem=128GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for discovery of CNV using WisecondorX
## Set hard coded paths and variables ##

modList=("Anaconda3/2024.06-1")
usage()
{
echo "# Script for discovery of CNV using WisecondorX.
# Requires: WisecondorX conda environment please see the readme for instructions on how to set this up.
#
# Usage sbatch --array 0-\"number of test samples -1\" -i sample.file.txt [ -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields: SampleName, SampleType , Sex.
# -c    OPTIONAL. Path to the configuration file.
# -r    REQUIRED. The name of your input reference file.
# -o    OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/ is used)
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
        -r )        shift
                    refFile=$1
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

if [ -z "${refFile}" ]; then # If no reference file name specified then epic failure.
    usage
    echo "## ERROR: The reference file was not specified this must be supplied."
    exit 1
fi
if [ ! -f "${refFile}" ]; then # If no reference file exists then epic failure.
    usage
    echo "## ERROR: The reference file ${refFile} was not found. Check the file exists and you do not have a typo in the name of the file supplied."
    exit 1
fi

readarray -t rSample <<< $(grep -v "^#" ${sampleFile} | grep test | cut -f1)
readarray -t Sex <<< $(grep -v "^#" ${sampleFile} | grep test | cut -f3)

if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir=${userDir}/variants/SV/wisecondorx
    echo "## INFO: Using ${outDir} as the output directory."
fi

if [ ! -d "${outDir}/cnvs" ]; then
    mkdir -p ${outDir}/cnvs
    echo "## INFO: Created the ${outDir}/cnvs directory."
fi

conda activate WisecondorX
wisecondorx predict \
--plot \
--add-plot-title  \
--bed \
--gender ${Sex[SLURM_ARRAY_TASK_ID]} \
${outDir}/wc/${Sample[SLURM_ARRAY_TASK_ID]}.merge.dedup.npz \
${refFile} \
${outDir}/cnvs/${Sample[SLURM_ARRAY_TASK_ID]}
conda deactivate
