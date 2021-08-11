#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=100:00:00     # Request runtime (up to 240 hours)
#$ -l h_vmem=10G     # Request RAM per core
#$ -m bea     # Status emails



module load python

python groupfile_reducing.py Galaxy12Filtersequencesbylengthondata11.fasta stability.contigs.groups renamed.groups




