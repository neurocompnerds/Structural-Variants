#!/bin/bash
#SBATCH -J clinSV
#SBATCH -o /hpcfs/users/%u/log/clinsv-slurm-%j.out

#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=18:00:00
#SBATCH --mem=144GB

# Notification Configuration 
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# This is the master script that coordinates job submission for the neurogenetics ClinSV structural variant pipeline.
## Set hard-coded paths and define functions ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
modList=("Singularity/3.10.5")

usage()
{
echo "# This script runs analysis of one or many Illumina whole genomes using ClinSV.
# Requires: BAM or CRAM files, Singularity and the ClinSV software. 
# Note: All output is written to /hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/test_run
#
# Usage sbatch $0 [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -o	OPTIONAL. /path/to/output. A default output directory will be used if this is not specified.
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 06/07/2021 
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
        -o )            shift
                        outDir=$1
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
if [ -z "${Config}" ]; then # If no config file specified use the default
    Config=${scriptDir}/configs/hs38DH.SV_ClinSV.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi

source ${Config}

if [ ! -f "${neuroDir}/clinsv/phoenix_resources.json" ]; then # Check if the config file exists
    cp "${defaultResourcesJSON}" "${neuroDir}/clinsv/phoenix_resources.json"
fi
if [ ! -d "${neuroDir}/clinsv/test_run" ]; then # Check if the output directory exists
    mkdir -p "${neuroDir}/clinsv/test_run"
fi
if [ -z "${outDir}" ]; then # If no output directory is specified use the default
    outDir=${neuroDir}/variants/SV/clinsv/clinsv_$(date +%Y%m%d_%H%M%S)
    echo "## INFO: Using the default output directory ${outDir}"
fi
if [ ! -d "${outDir}" ]; then # Check if the output directory exists
    mkdir -p "${outDir}"
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load ${mod}
done

## Launch ClinSV ##
singularity exec --bind ${refPath}:/app/ref-data/refdata-b38,${neuroDir}/clinsv/test_run:/app/project_folder,${neuroDir}/clinsv:/app/input \
${progDir}/${progName} /app/clinsv/bin/clinsv \
-p /app/project_folder/ \
-ref /app/ref-data/refdata-b38 \
-i "/app/input/*.bam" \
-j /app/input/phoenix_resources.json \
-w 

# Clean up and reset for the next run
if [ -f "${neuroDir}/clinsv/test_run/results/SV-CNV.RARE_PASS_GENE.xlsx"]; then # Proxy test for run completion
    rm ${neuroDir}/clinsv/*.bam ${neuroDir}/clinsv/*.bai # Remove the input BAM files to save space
else
    echo "## ERROR: ClinSV did not produce the expected output file. Looks like something went wrong you may need to check the logs for errors."
    rm "${neuroDir}/clinsv/clinsv.lock" # Remove the lock file to allow other runs to proceed
    exit 1
fi
mv ${neuroDir}/clinsv/test_run/* "${outDir}/" # Move the output files to the output directory
echo "## INFO: ClinSV run completed. Results are in ${outDir}"

if [ -f "${neuroDir}/clinsv/clinsv.lock" ]; then # Check if the lock file exists
    rm "${neuroDir}/clinsv/clinsv.lock" # Remove the lock file to allow other runs to proceed
fi
