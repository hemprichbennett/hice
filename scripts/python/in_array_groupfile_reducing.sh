#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=34:00:00     # Request runtime (up to 240 hours)
#$ -l h_vmem=5G     # Request RAM per core
#$ -t 1-500





file=$(ls x* | sed -n ${SGE_TASK_ID}p)

echo working with $file

module load python

echo loaded python

python array_groupfile_reducing.py $file stability.contigs.groups 




