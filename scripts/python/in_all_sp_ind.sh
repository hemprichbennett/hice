#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=70:0:0     # Request runtime (up to 240 hours)
#$ -l h_vmem=75G     # Request RAM per core
#$ -m bea     # Status emails



module load R



Rscript /data/scratch/btw863/bat-diet/scripts/r/all_species_network_generation.R
