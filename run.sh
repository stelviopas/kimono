#!/bin/bash

#SBATCH -o miss_data_imput.out
#SBATCH -e miss_data_imput.err
#SBATCH -J miss_data_imput
#SBATCH -p cpu_p
#SBATCH -c 4 # default values is 2
#SBATCH --mem=16G # total memory in MB
#SBATCH -t 10:00:00 # format dd-hh:mm:ss
#SBATCH --nice=10000  # adjusts scheduling priority

echo "Starting batch job at `date`"
Rscript ~/kimono/R/impute.R

