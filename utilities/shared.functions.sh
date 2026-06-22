#!/bin/bash

test_genome_build() {
case "${Build}" in
    "hs38DH" | "GRCh38_full_analysis_set" )    genomeType="has_alt_contigs"
                ;;
    "hs37d5" | "ucsc.hg19" )    genomeType="no_alt_contigs"
                ;;
    "CHM13v2" )    genomeType="no_alt_contigs"
                ;;
    * )         echo "## ERROR: The genome build ${Build} was not recognized: Available options are: hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 and CHM13v2"
                echo "You can add new genomes by editing the test_genome_build function in the file configs/BWA-GATKHC.environment.cfg"
                echo "You should also create a config file and the appropriate GATK reference files for your new genome similar to those in the configs/ directory"
                exit 1
                ;;
esac
}

select_genome_build()
{
case "${genomeSize}" in
    3099922541 )    buildID="GRCh38"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                    ;;
    3217346917 )    buildID="hs38DH"
                    genomeBuild="$refDir/hs38DH.fa"
                    ;;
    3137454505 )    buildID="hs37d5"
                    genomeBuild="$refDir/hs37d5.fa.gz"
                    ;;
    2730871774 )    buildID="GRCm38"   
                    genomeBuild="$refDir/GRCm38_68.fa"
                    ;;
    3117463893 )    buildID="CHM13v2"
                    genomeBuild="$refDir/T2T_CHM13v2.0.ucsc.ebv.fa.gz"
                    ;;
    3137161264 )    buildID="hg19"
                    genomeBuild="$refDir/ucsc.hg19.fasta"
                    ;;
    3105715063 )    buildID="GRCh38.hs38d1.no_alt"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3215250450 )    buildID="GRCh38.hs38d1.full"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3099750718 )    buildID="GRCh38"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                    ;;
    3031042417 )    buildID="GRCh38.blacklist"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
                    ;;
    * )         echo "## ERROR: Genome length $genomeSize for ${bamFile[SLURM_ARRAY_TASK_ID]} was not matched, you may need to specify the genome build directly using the -g flag."
                exit 1
                ;;
esac
echo "## INFO: The CRAM file ${bamFile[SLURM_ARRAY_TASK_ID]} was likely mapped to $buildID corresponding to the refseq $genomeBuild."
}
