#!/bin/bash

#SBATCH -J retroseq
#SBATCH -o /fast/users/%u/retroseq-slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
usage()
{
echo "# Script for identifying retrotransposed elements in an aligned genome using RetroSeq pipeline
# A genome usually requires 16 cores and time 1 day, mem 8GB
#
# Usage sbatch --array 0-(# of bam files minus 1) $0 -i /path/to/input [-a percent-align-identity -o /path/to/output -g /path/to/genome.fa -r /path/to/repeatsDir ] | [ - h | --help ]
#
# Options
# -i    REQUIRED. Location of your bam files to run retroseq on
# -a	OPTIONAL. Equivalent to RetroSeq -id option (NOTE: program default is 90 but default for this script is 80)
# -o    OPTIONAL. Path to where you want to find your file output (if not specified $FASTDIR/retroseq is used)
# -g    OPTIONAL. Path to your reference genome (if not specified /data/neurogenetics/RefSeq/ucsc.hg19.fasta is used)
# -r    OPTIONAL. Path to where repeats fasta files are (if not specified /data/neurogenetics/RefSeq/repeats is used)
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
#
# Original: Atma Ivancevic, 20/08/2018
# Modified: (Date; Name; Description)
# 27/04/2019; Mark Corbett <mark dot corbett at adelaide.edu.au>; Added usage and script flags; centralised perl script and default locations of repeats directory
#
"
}

## Set default program location
RETROSEQEXE=/data/neurogenetics/executables/RetroSeq/bin/retroseq.pl

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -i )                    shift
                                        INDIR=$1
                                        ;;
                -a )                    shift
                                        ALIGNID=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                -g )                    shift
                                        GENOME=$1
                                        ;;
                -r )                    shift
                                        REPEATSDIR=$1
                                        ;;
                -h | --help )           usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

if [ -z "$INDIR" ]; then # If no input directory specified then do not proceed
        usage
        echo "#ERROR: You need to tell me where to find some bam files to run the script on are."
        exit 1
fi
if [ -z "$ALIGNID" ]; then # If no alignment minimum match percentage specified use the default
        ALIGNID=80
fi
if [ -z "$OUTDIR" ]; then # If no output directory then use the default
        OUTDIR=$FASTDIR/retroseq
fi
if [ -z "$GENOME" ]; then # If no genome specified then use the default
        GENOME=/data/neurogenetics/RefSeq/ucsc.hg19.fasta
fi
if [ -z "$REPEATSDIR" ]; then # If no repeat directory then use the default
        REPEATSDIR=/data/neurogenetics/RefSeq/repeats
fi

if [ ! -d $OUTDIR ]; then
        mkdir -p $OUTDIR
fi
if [ ! -d $OUTDIR/TEdiscovery ]; then
        mkdir -p $OUTDIR/TEdiscovery
        mkdir -p $OUTDIR/TEcalling
fi

BUILD=$(basename $GENOME)

# define query bam files
mapfile -t QUERIES < <(find $INDIR/*.bam | xargs -n1 basename)

# Do a little summary for the log file
echo "
# Running Retroseq with the following parameters
# File: ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} 
# Outputs can be found here: $OUTDIR
# Minimum match percentage: $ALIGNID %
# Genome: $BUILD (If you have problems with the output make sure this reference matches the BAM file)
# Repeat library: $REPEATSDIR
"
# load modules
module load Exonerate/2.2.0-foss-2016uofa
module load BEDTools/2.25.0-GCC-5.3.0-binutils-2.25
module load SAMtools/0.1.19-GCC-5.3.0-binutils-2.25

# run pipeline

### discovery phase ###
echo "discovering..."
$RETROSEQEXE -discover \
-bam ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
-output ${OUTDIR}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab \
-eref ${REPEATSDIR}/eref_types.tab \
-align -id ${ALIGNID} > ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log 2>&1
echo "done"

### filtering ###
# filter out candidates near contig starts (likely contamination)
echo "filtering contigs..."
awk -F"\t" '{if ($2>1000) print}' ${OUTDIR}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab \
> ${OUTDIR}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab.filtered
echo "done"

### calling phase ###
echo "calling..."
$RETROSEQEXE -call \
-bam ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
-input ${OUTDIR}/TEdiscovery/${QUERIES[$SLURM_ARRAY_TASK_ID]}.candidates.tab.filtered \
-ref ${GENOME} \
-output ${OUTDIR}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf >> ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log 2>&1
echo "done"

module load HTSlib/1.9-foss-2016b
bgzip ${OUTDIR}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf && tabix ${OUTDIR}/TEcalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf.gz

# Collect the slurm log file
mv $FASTDIR/retroseq-slurm-$SLURM_JOB_ID.out $OUTDIR/
