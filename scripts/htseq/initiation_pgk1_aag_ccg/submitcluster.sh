#!/bin/bash -l

#SBATCH --mem=32000
#SBATCH -t 00:30:00  

# runs this specific script on cluster
python calculate_unique_insert_counts.py $1
