#!/bin/bash
# This script coordinates submission of sequential jobs for WISECONDORX.
## Hard coded paths for your system should be set in configs/BWA-GATKHC.environment.cfg ##

usage()
{
echo "# This is the master script that coordinates job submission for analysis of low coverage genome sequencing data using the WISECONDORX pipeline.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# The scripts deliver a sorted and deduplicated BAM file, a WGS metrics report, a pipeline log and WISCONDORX reports.
# Requires: BWA-MEM, samtools, sambamba, GATK / Picard, Java, WisecondorX conda environment.
#
# Usage $0 -i input.file.txt [ -c /path/to/Config.cfg -o /path/to/output ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. A path to a tab delimited text file that has the required fields SampleName, /path/to/R1.fastq.gz, /path/to/R2.fastq.gz, sample_type (test or control). Plus optional fields for library, platform and platform unit.
# -c    OPTIONAL. Path to the configuration file.
# -o    OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
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

nTasks=$(($(grep -vc "^#" ${inputFile})-1)) # Count the number of non-comment lines in the input file to get the number of read pairs minus one.

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

if [ -z "${outDir}" ]; then # If no output directory then set and create a default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi
if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
fi
for dir in qc align wc cnvs; do
    if [ ! -d "${outDir}/${dir}" ]; then
        mkdir -p ${outDir}/${dir}
    fi
done

## Make a file that just has the sample name and sample type columns
sampleFile=${outDir}/$(date +"[%Y%m%d]").wisecondorx.samples.txt
cut -f1,4,5 ${inputFile} | grep -v "^#" | sort | uniq > ${sampleFile}
nSamples=$(wc -l ${sampleFile}) # Count the number of unique sample names in the input file to get the number of samples.
nSampleTasks=$(($nSamples-1)) # Get the number of samples minus one for array job submission.
nTestSamples=$(grep test ${sampleFile} | wc -l)
nTestSampleTasks=$((${nTestSamples}-1))
refFile=${outDir}/$(date +"[%Y%m%d]").wisecondorx.reference.npz


## Launch the job chain ##
# A wiseman, a con artist and two dorks walk into a bar.  The wiseman says "Myrrh! I should have seen that coming!" The con artist says "I'll sell you some snake oil to make you feel better" one the two dorks said "Ahhh my nuts!" and the other one laughed.
mapJob=`sbatch --array=0-${nTasks} --export=ALL mapSort_alt_aware.sh -i ${inputFile} -c ${enviroCfg} -o ${outDir}` 
mapJob=$(echo ${mapJob} | cut -d" " -f4)
mergeJob=`sbatch --array=0-${nSampleTasks} --export=ALL mergeBam.sh --dependency=afterok:${mapJob} -i ${sampleFile} -c ${enviroCfg} -o ${outDir}`
mergeJob=$(echo ${mergeJob} | cut -d" " -f4)
mosdepthJob=`sbatch --array=0-${nSampleTasks} --export=ALL --dependency=afterok:${mergeJob} mosdepth.sh -i ${sampleFile} -c ${enviroCfg} -o ${outDir}`
mosdepthJob=$(echo ${mosdepthJob} | cut -d" " -f4)
convertJob=`sbatch --array=0-${nSampleTasks} --export=ALL --dependency=afterok:${mergeJob} wisecondorx.convert.sh -i ${sampleFile} -c ${enviroCfg} -o ${outDir}`
convertJob=$(echo ${convertJob} | cut -d" " -f4)
newRefJob=`sbatch --export=ALL --dependency=afterok:${convertJob} wisecondorx.create.ref.sh -i ${sampleFile} -r ${refFile} -c ${enviroCfg} -o ${outDir}`
newRefJob=$(echo ${newRefJob} | cut -d" " -f4)
predictJob=`sbatch --array=0-${nTestSampleTasks} --export=ALL --dependency=afterok:${newRefJob} wisecondorx.predict.sh -i ${sampleFile} -r ${refFile} -c ${enviroCfg} -o ${outDir}`
predictJob=$(echo ${predictJob} | cut -d" " -f4) # Placeholder until further development on the QC outputs
