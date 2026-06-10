#!/bin/bash

test_genome_build() {
case "${BUILD}" in
    "hs38DH" | "GRCh38_full_analysis_set" )    genomeType="has_alt_contigs"
                ;;
    "hs37d5" | "ucsc.hg19" )    genomeType="no_alt_contigs"
                ;;
    "CHM13v2" )    genomeType="no_alt_contigs"
                ;;
    * )         echo "## ERROR: The genome build ${BUILD} was not recognized: Available options are: hs38DH, GRCh38_full_analysis_set, hs37d5, ucsc.hg19 and CHM13v2"
                echo "You can add new genomes by editing the test_genome_build function in the file configs/BWA-GATKHC.environment.cfg"
                echo "You should also create a config file and the appropriate GATK reference files for your new genome similar to those in the configs/ directory"
                exit 1
                ;;
esac
}
