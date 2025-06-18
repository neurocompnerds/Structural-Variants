# Structural-Variants

## Introduction
Our collection of scripts used to provide inputs to tools that detect structural variants and/or novel sequence insertions in the human genome. The scripts are designed specifically for use on the University of Adelaide HPC and for our Neurogenetics research program. Others may be able to use our general framework but it will need tweaking especially of the config files to work in your hands. I'd suggest checking out nf-core which probably has a pipeline that does the same but with better reporting and QC.

This repo doesn't get as much love as it should so even if you are from the Neuro team you may find some scripts do not work on the current set up of the Adelaide University HPC.

Scripts / pipelines that *should* work:

GRIDSS

Parliament2

MEI/retroseq

clinsv

Things that definitely *don't* work:

delly (use Parliament2 instead).

smoove

## Usage notes
Some of these pipelines are described in our Google Doc protocols.

https://docs.google.com/document/d/1TbvCCqPcnCx9Xl9MZPZgQquyuG6FG43shAobZbsFB8I

Most scripts default to hs38DH build which is also the default build used for our map-n-call (BWA/GATK) pipeline https://github.com/speleonut/map-n-call.  If you want to use these tools on a different genome build then you will need to make a config file for it (see 'configs/' directory for examples) *and* you need to make sure that all of the neccesary reference files are available for it.

Most scripts will have help built-in and these can usually be accessed with `-h` or `--help` flags when executing the script on the command line.  

GRIDSS is the most complicated pipeline which has multiple stages associated with it.  You can execute the entire pipeline from the command line using the `GRIDSS/gridss_launcher.sh` script which will set up all of your SLURM submissions.

Parliament2 is not kind to your available disk space.  It is limited due to hard coded paths and file names that it uses and generalised locations for dumping of temporary files while it runs different programs.  To get it to run in parallel without each sample overwriting files from another sample this script needs to copy the genome reference to a folder for *every* BAM/CRAM file you are running.  Our script makes a copy of the input BAM/CRAM file because Parliament2 also does not respect input data (it will move it, rename it and delete the original). The script makes further replications of the BAM when it splits everything by chromosome.  Future incarnations of the script will do a better job of cleaning up after the process has run.

The ClinSV pipeline can analyse multiple genomes at once. It is strongly advised to keep it to just a few genomes at a time due to the high requirement for disk space. In future we will edit the pipeline to submit either low or high resource requests depending on the number of genomes submitted to aid in getting jobs to start quickly.  Note that there can be only one ClinSV run at a time as all output must be run from `/hpcfs/groups/phoenix-hpc-neurogenetics/clinsv/`. When the pipeline is running, a "lock file" will prevent you from starting a second ClinSV job.  If the pipeline fails to complete then the lock file may not be properly cleared however it can be easily removed manually (the script will tell you how to do this).  You can check if the pipeline is running by searching for the job in the queue e.g. `squeue | grep clinSV`.  When you have finished your 

