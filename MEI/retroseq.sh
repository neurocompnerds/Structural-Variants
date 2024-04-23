#!/bin/bash

#SBATCH -J retroseq
#SBATCH -o /hpcfs/users/%u/retroseq-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
neuroDir=/hpcfs/groups/phoenix-hpc-neurogenetics
retroseqScript=$neuroDir/executables/RetroSeq/bin/retroseq.pl
modList=("arch/haswell" "Exonerate/2.2.0-foss-2016uofa" "BEDTools/2.25.0-GCC-5.3.0-binutils-2.25" "SAMtools/0.1.19-GCC-5.3.0-binutils-2.25")

usage()
{
echo "# Script for identifying retrotransposed elements in an aligned Genome using RetroSeq pipeline
# A Genome usually requires 16 cores and time 1 day, mem 8GB
#
# Usage sbatch --array 0-(# of bam files minus 1) $0 -i /path/to/input [-a percent-align-identity -o /path/to/output -g /path/to/Genome.fa -r /path/to/repeatsDir ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. Location of your bam files to run retroseq on
# -a	OPTIONAL. Equivalent to RetroSeq -id option (NOTE: program default is 90 but default for this script is 80)
# -o    OPTIONAL. Path to where you want to find your file output (if not specified /hpcfs/users/${USER}/retroseq is used)
# -g    OPTIONAL. Path to your reference Genome (if not specified $neuroDir/RefSeq/ucsc.hg19.fasta is used)
# -r    OPTIONAL. Path to where repeats fasta files are (if not specified $neuroDir/RefSeq/repeats is used)
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
#
# Original: Atma Ivancevic, 20/08/2018
# Modified: (Date; Name; Description)
# 27/04/2019; Mark Corbett <mark dot corbett at adelaide.edu.au>; Added usage and script flags; centralised perl script and default locations of repeats directory
#
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -i )                    shift
                                        inputDir=$1
                                        ;;
                -a )                    shift
                                        pctAlignId=$1
                                        ;;
                -o )                    shift
                                        outputDir=$1
                                        ;;
                -g )                    shift
                                        Genome=$1
                                        ;;
                -r )                    shift
                                        repeatsDir=$1
                                        ;;
                -h | --help )           usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

if [ -z "$inputDir" ]; then # If no input directory specified then do not proceed
    usage
    echo "## ERROR: You need to tell me where to find some bam files to run the script on are."
    exit 1
fi
if [ -z "$pctAlignId" ]; then # If no alignment minimum match percentage specified use the default
    pctAlignId=80
fi
if [ -z "$outputDir" ]; then # If no output directory then use the default
    outputDir=/hpcfs/users/${USER}/retroseq
    echo "## INFO: No output directory was specified so the default $outputDir will be used."
fi
if [ -z "$Genome" ]; then # If no Genome specified then use the default
    Genome=$neuroDir/RefSeq/ucsc.hg19.fasta
fi
if [ -z "$repeatsDir" ]; then # If no repeat directory then use the default
        repeatsDir=$neuroDir/RefSeq/repeats
fi

if [ ! -d $outputDir ]; then
        mkdir -p $outputDir
fi
if [ ! -d $outputDir/TEdiscovery ]; then
        mkdir -p $outputDir/TEdiscovery
        mkdir -p $outputDir/TEcalling
fi

BUILD=$(basename $Genome)

# define query bam files
mapfile -t QUERIES < <(find $inputDir/*.bam | xargs -n1 basename)

# Do a little summary for the log file
echo "
# Running Retroseq with the following parameters
# File: ${inputDir}/${QUERIES[$SLURM_ARRAY_TASK_ID]} 
# Outputs can be found here: $outputDir
# Minimum match percentage: $pctAlignId %
# Genome: $BUILD (If you have problems with the output make sure this reference matches the BAM file)
# Repeat library: $repeatsDir
"
# load modules

# run pipeline

### discovery phase ###
echo "discovering..."

retroseqScript -discover \
-bam ${inputDir}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
-output ${outputDir}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab \
-eref ${repeatsDir}/eref_types.tab \
-align -id ${pctAlignId} > ${outputDir}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log 2>&1
echo "done"

### filtering ###
# filter out candidates near contig starts (likely contamination)
echo "filtering contigs..."
awk -F"\t" '{if ($2>1000) print}' ${outputDir}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab \
> ${outputDir}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab.filtered
echo "done"

### calling phase ###
echo "calling..."

retroseqScript -call \
-bam ${inputDir}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
-input ${outputDir}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab.filtered \
-ref ${Genome} \
-output ${outputDir}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf >> ${outputDir}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log 2>&1
echo "done"

module load HTSlib/1.9-foss-2016b
bgzip ${outputDir}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf && tabix ${outputDir}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf.gz

# Collect the slurm log file
mv /hpcfs/users/${USER}/retroseq-slurm-$SLURM_JOB_ID.out $outputDir/
