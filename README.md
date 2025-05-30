# Structural-Variants

## Introduction
Our collection of scripts used to provide inputs to tools that detect structural variants and/or novel sequence insertions in the human genome. The scripts are designed specifically for use on the University of Adelaide HPC and for our Neurogenetics research program.  If you aren't using that then you may be able to use our general framework but it will need tweaking especially of the config files to work in your hands. I'd suggest checking out nf-core which probably has a pipeline that does the same but with better reporting and QC.

This repo doesn't get as much love as it should so even if you are from the Neuro team you may find some scripts do not work on the current set up of the Adelaide University HPC.

Scripts / pipelines that *should* work:
GRIDSS
Parliament2
MEI/retroseq

Things that definitely *don't* work:
clinsv (we wish it would).
delly (use Parliament2 instead).
smoove

## Usage notes
Most scripts default to hs38DH build which is the default build also used for our map-n-call (BWA/GATK) pipeline.  If you want to use these tools on a different genome build then you need to make a config file for it (see 'configs/' directory for examples) *and* you need to make sure that all of the neccesary reference files are available for it.

Most scripts will have help built-in and these can usually be accessed with '-h' or '--help' flags when executing the script on the command line.  

GRIDSS is the most complicated pipeline which has multiple stages associated with it.  You can execute the entire pipeline from the command line using the 'GRIDSS/gridss_launcher.sh' script which will set up all of your SLURM submissions.

