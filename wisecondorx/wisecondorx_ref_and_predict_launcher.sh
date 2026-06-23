#!/bin/bash
# This script coordinates submission of sequential jobs for WISECONDORX.
## Hard coded paths for your system should be set in configs/BWA-GATKHC.environment.cfg ##

usage()
{
echo "# This is the short-cut script that coordinates job submission for analysis of low coverage genome sequencing data using the WISECONDORX pipeline.
# The script takes a tab-delimited input file and submits jobs for reference creation and prediction.
# The scripts deliver WISCONDORX reports.
# Requires: WisecondorX conda environment.
#
# Usage $0 -i input.file.txt [ -c /path/to/Config.cfg -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields SampleName, sample_type (test or control) and Sex. You can use the input file from the full workflow if you like. 
# -c    OPTIONAL. Path to the configuration file.
# -o    REQUIRED. Path to where your outputs from running the full pipeline are located (you need the .npz files for each sample to run this script).
# -h or --help    Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 03/06/2026 
# Modified: (Date; Name; Description)
# See https://github.com/speleonut/map-n-call/activity for script edit history
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
        -h | --help )   usage
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
if [ -z "${outDir}" ]; then # If no output directory then set and create a default directory
    usage
    echo "## ERROR: You need to tell me where your samples are located.
    # -o    REQUIRED. Path to where your outputs from running the full pipeline are located (you need the .npz files for each sample to run this script)."
    exit 1
fi
if [ ! -d "${outDir}" ]; then
    echo "## ERROR: The output directory ${outDir} was not found."
    exit 1
fi

if [ -z "${enviroCfg}" ]; then # If there is no config file then use the default.
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,wisecondorx,configs,g')"
    enviroCfg="${configDir}/hs38DH.SV_WISECONDORX.phoenix.cfg"
    echo "## INFO: the default config file ${enviroCfg} has been selected."
fi

source ${enviroCfg}

if [ -z "${Build}" ]; then # If Build is not in the config then let the user know and exit
    echo "## ERROR: Genome Build not specified in the config file.  Please check your config file and make sure the Build variable is set to one of the following: hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 or CHM13v2."
    exit 1
fi

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi
if [ ! -d "${baseTmpDir}" ]; then
    mkdir -p ${baseTmpDir}
fi

## Check the input file and if needed make a file that just has the sample name, sample type and sex columns
dateString=$(date +"%Y%m%d_%H%M%S")
sampleFile=${outDir}/${dateString}.wisecondorx.samples.txt
testInput=$(grep -v "^#" ${inputFile} | cut -f4 | head -n1)
if [[ "${testInput}" != "test" && "${testInput}" != "control" ]]; then
    echo "## INFO: Testing your input file suggests it is likely to be an input file from the full pipeline."
    grep -A1 "#Sample" ${inputFile}
    echo "## INFO: Creating a sample file ${sampleFile} for this workflow."
    cut -f1,4,5 ${inputFile} | grep -v "^#" | sort | uniq > ${sampleFile}
else
    sampleFile=${inputFile}
fi
nSamples=$(grep -cv "^#" ${sampleFile}) # Count the number of unique sample names in the input file to get the number of samples.
nSampleTasks=$(($nSamples-1)) # Get the number of samples minus one for array job submission.
nTestSamples=$(grep test ${sampleFile} | wc -l)
nTestSampleTasks=$((${nTestSamples}-1))
refFile=${outDir}/${dateString}.wisecondorx.reference.npz

## Launch the job chain ##
newRefJob=`sbatch --export=ALL ${scriptDir}/wisecondorx/wisecondorx.create.ref.sh -i ${sampleFile} -r ${refFile} -c ${enviroCfg} -o ${outDir}`
newRefJob=$(echo ${newRefJob} | cut -d" " -f4)
predictJob=`sbatch --array=0-${nTestSampleTasks} --export=ALL --dependency=afterok:${newRefJob} ${scriptDir}/wisecondorx/wisecondorx.predict.sh -i ${sampleFile} -r ${refFile} -c ${enviroCfg} -o ${outDir}`
predictJob=$(echo ${predictJob} | cut -d" " -f4) # Placeholder until further development on the QC outputs
