# WISECONDORX 
This software is for determining copy number variants from low coverage short-read genome sequencing.  To use it on Phoenix you will need to set up the conda environment.
https://github.com/CenterForMedicalGeneticsGhent/WisecondorX

See the HPC wiki for basic set up of conda if you don't already have this (you'll need to login with your Adelaide University credentials to view this):
https://wiki.adelaide.edu.au/hpc/Anaconda  If it is a while since you set this up then you might want to clear out all of the old conda paths from your `.bashrc` and other `.*rc` files and re run `conda init`.  Use `conda config --show` to check everything is pointing to Anaconda3/2024.06-1. Don't forget to log out and back in so that the changes are registered (I guess you could also `source` the files if needed).

Set up your WisecondorX environment (you'll need to use the same case as indicated below to get this pipeline to work for you).

```bash
ml Anaconda3/2024.06-1
conda create --name WisecondorX -c conda-forge -c bioconda wisecondorx
```

You should be good to go now.
