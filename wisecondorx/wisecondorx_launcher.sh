#!/bin/bash
# This script coordinates submission of sequential jobs for WISECONDORX.
## Hard coded paths for your system should be set in configs/BWA-GATKHC.environment.cfg ##
whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
configDir="$(echo ${whereAmI} | sed -e 's,wisecondorx,configs,g')"
enviroCfg="${configDir}/hs38DH.SV_WISECONDORX.phoenix.cfg"
source ${enviroCfg}

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

if [ ! -d "${baseTmpDir}" ]; then
	mkdir -p ${baseTmpDir}
fi

usage()
{
echo "# This is the master script that coordinates job submission for primarily Illumina genome sequencing alignments but will work for exomes too.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# The scripts deliver an indel realigned BAM file, a WGS metrics report, a pipeline log and a gzipped gVCF file.
# Requires: BWA-MEM, samtools, sambamba, GATK / Picard, Java, WisecondorX conda environment.
#
# Usage $0 -i input.file.txt [ -g genomeBuild-o /path/to/output ] | [ - h | --help ]
#
# Options
# -i	REQUIRED. A path to a tab delimited text file that has the required fields SampleName, /path/to/R1.fastq.gz, /path/to/R2.fastq.gz, sample_type (test or control). Plus optional fields for library, platform and platform unit.
# -g	OPTIONAL. A code for a genome build, either hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 or CHM13v2.  If not specified then the default is hs38DH.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
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
		-i )			shift
					inputFile=$1
					;;
		-g )			shift
					genomeBuild=$1
					;;
		-o )			shift
					outDir=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
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
nSamples=$(cut -f1 ${inputFile} | grep -v "^#" | sort | uniq | wc -l) # Count the number of unique sample names in the input file to get the number of samples.

if [ -z "${outDir}" ]; then # If no output directory then set and create a default directory
	outDir="${userDir}/variants/SV/wisecondorx/"
	echo "## INFO: Using ${outDir} as the output directory"
fi
if [ ! -d "${outDir}" ]; then
	mkdir -p ${outDir}
fi
if [ -z "${genomeBuild}" ]; then # If no genome build specified, set to default
    genomeBuild=hs38DH
    echo "## INFO: No genome build specified, using ${genomeBuild} as the default"
fi

## Launch the job chain ##
# A wiseman, a con artist and two dorks walk into a bar.  The wiseman says "Myrrh! I should have seen that coming!" The con artist says "I'll sell you some snake oil to make you feel better" one the two dorks said "Ahhh my nuts!" and the other one laughed.
