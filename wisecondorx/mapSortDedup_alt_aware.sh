#!/bin/bash

#SBATCH -J MapSortDup
#SBATCH -o /hpcfs/users/%u/log/mapSortDedup-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --time=06:00:00
#SBATCH --mem=148GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for mapping Illumina pair-end sequence data
## Set hard coded paths and variables ##
source ${enviroCfg}

modList=("BWA/0.7.17-GCCcore-11.2.0" "HTSlib/1.17-GCC-11.2.0" "SAMtools/1.17-GCC-11.2.0")
usage()
{
echo "# Script for mapping Illumina pair-end sequence data
# Requires: BWA 0.7.x, samtools, sambamba
# This script assumes your sequence files are gzipped
#
# Usage $0 -i input.file.txt [ -g genomeBuild-o /path/to/output ] | [ - h | --help ]
#
# Options
# -i	REQUIRED. A path to a tab delimited text file that has the required fields SampleName, /path/to/R1.fastq.gz, /path/to/R2.fastq.gz. Plus optional fields for library, platform and platform unit.  Please use the provided template or unpredictable results and/or embarrassing failure is assured.
# -g	OPTIONAL. A code for a genome build, either hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 or CHM13v2.  If not specified then the default is hs38DH.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory ${userDir}/variants/SV/wisecondorx/\${Sample} is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 03/06/2026; Mark Corbett; Fork from https://github.com/speleonut/map-n-call/blob/master/GATK4/mapSortDedupMarkIndels_alt_aware_phoenix.sh and converted to an array-capable script.
#
"
}

## Set Variables ##
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

if [ -z "${outDir}" ]; then # If no output directory then use the default directory
    outDir="${userDir}/variants/SV/wisecondorx/"
    echo "## INFO: Using ${outDir} as the output directory"
fi

if [ -z "${genomeBuild}" ]; then # If no genome build specified, set to default
    genomeBuild=hs38DH
    echo "## INFO: No genome build specified, using ${genomeBuild} as the default"
fi

test_genome_build # This function is defined in the configs/ file and sets the variable genomeType to either "has_alt_contigs" or "no_alt_contigs" depending on the genome build specified.

## Parse the input file ##
# This assumes you use the provided template and keep all of the columns in the right order.

Sample=$((grep -v "^#" ${inputFile} | cut -f1)) 
seqFile1=$((grep -v "^#" ${inputFile} | cut -f2)) 
seqFile2=$((grep -v "^#" ${inputFile} | cut -f3)) 
LB=$(grep -v "^#" ${inputFile} | cut -f4) # Library information (optional)
PL=$(grep -v "^#" ${inputFile} | cut -f5) # Platform information (optional)
PU=$(grep -v "^#" ${inputFile} | cut -f6) # Platform unit information (optional)
ID=$(zcat ${seqFile1[SLURM_ARRAY_TASK_ID]} | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').$outPrefix # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq

## Create essential directories ##

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
    echo "## INFO: ${outDir} has been created as the output directory"
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
# Map reads to genome using BWA-MEM
# do NOT use the -M option when aligning to GRCh38 as alignment is alt-aware
# -K flag asks bwa-mem to load a fixed number of bases into RAM so enables reproducibility
 
cd $tmpDir
bwa mem -K 100000000 -t 24 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:$Sample" $BWAINDEXPATH/$BWAINDEX $seqFile1 $seqFile2 |\
samtools view -bT $GATKREFPATH/$BUILD/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@24 -T$Sample -o $Sample.samsort.bwa.$BUILD.bam -

# Mark duplicates
$sambambaProg markdup -t 24 -l 5 --tmpdir=$tmpDir --overflow-list-size 1000000 --hash-table-size 1000000 $Sample.samsort.bwa.$BUILD.bam $outDir/$Sample.marked.sort.bwa.$BUILD.bam
if [ -f "$outDir/$Sample.marked.sort.bwa.$BUILD.bam" ]; then
    rm $Sample.samsort.bwa.$BUILD.bam
else
	echo "## ERROR: Duplicate marking or earlier stage failed!"
	exit 1
fi

echo "# Flagstats" > $outDir/$Sample.Stat_Summary.txt
samtools flagstat $outDir/$Sample.marked.sort.bwa.$BUILD.bam >> $outDir/$Sample.Stat_Summary.txt

if [ -f "$tmpDir/${outPrefix}.cat_R2.fastq.gz" ]; then
    rm $tmpDir/${outPrefix}.cat_R1.fastq.gz $tmpDir/${outPrefix}.cat_R2.fastq.gz
fi