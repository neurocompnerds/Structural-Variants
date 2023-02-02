#!/bin/bash
#SBATCH -J assembleGRIDSS
#SBATCH -o /hpcfs/users/%u/log/assembleGRIDSS-slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=1-00:00:00
#SBATCH --mem=32GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=${USER}@adelaide.edu.au

# A script to preprocess bam files for GRIDSS
## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Structural-Variants"
customModDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/easybuild/modules/all"
modList=("arch/skylake" "SAMtools/1.12" "Java/1.8.0_191" "BWA/0.7.17" "R/4.0.3")
RLibDir="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/R/4.0.3/RLibs"
threads=8
assembly_jobs=32

usage()
{
echo "# This preprocesses BAM files for GRIDSS.  This script can be directly submitted to the scheduler but it is best to let gridss_launcher.sh handle this.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# Requires: An aligned BAM file (or files), BWA, samtools, Java, R 
#
# Usage sbatch $0 -p file_prefix -i /path/to/input/bam-file [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX*.bam
# -i	REQUIRED. Path to the BAM file
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/users/${USER}/GRIDSS/\${outPrefix} is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 06/07/2021 
# Modified: (Date; Name; Description)
# 
# 
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -c )            shift
                        Config=$1
                        ;;
        -p )            shift
                        outPrefix=$1
                        ;;
        -i )            shift
                        inputDir=$1
                        ;;
        -o )            shift
                        workDir=$1
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
if [ -z "$Config" ]; then # If no config file specified use the default
    Config=$scriptDir/configs/hs38DH.SV_GRIDSS.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "${outPrefix}" ]; then # If no file prefix specified then do not proceed
    usage
    echo "## ERROR: You need to specify a file prefix (PREFIX) referring to your BAM files eg. PREFIX*.bam."
    exit 1
fi
if [ -z "$inputDir" ]; then # If path to bam file not specified then do not proceed
	usage
	echo "## ERROR: You need to specify the path to your bam files"
	exit 1
fi
# Locate the bam
bamFile=$(find $inputDir/*.bam | grep $outPrefix)
if [ ! -f "$bamFile" ]; then # Check if the BAM is actually CRAM
    bamFile=$(find $inputDir/*.cram | grep $outPrefix)
    if [ ! -f "$bamFile" ]; then
        echo "## ERROR: BAM or CRAM file not found in $inputDir"
        exit 1
    fi
fi 
if [ -z "$workDir" ]; then # If no output directory then set and create a default directory
	workDir=/hpcfs/users/${USER}/GRIDSS/$outPrefix
	echo "## INFO: Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
    mkdir -p $workDir
fi
tmpDir=$baseTmpDir/$outPrefix
if [ ! -d "$tmpDir" ]; then
    mkdir -p $tmpDir
fi

## Load modules ##
module use $customModDir
for mod in "${modList[@]}"; do
    module load $mod
done
export R_LIBS_USER=${RLibDir}

## Create a job array containing assembly jobs

$gridss_cmd_common \
-t $threads \
-o $workDir/$outPrefix.sv.vcf.gz \
-a $workDir/$outPrefix.asm.bam \
-w $tmpDir \
-s assemble \
--jobindex $SLURM_ARRAY_TASK_ID \
--jobnodes $assembly_jobs \
$bamFile >> $workDir/$outPrefix.gridss.log 2>&1
