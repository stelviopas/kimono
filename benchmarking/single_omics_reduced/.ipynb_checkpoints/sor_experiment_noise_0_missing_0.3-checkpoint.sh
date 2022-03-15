#!/bin/bash

#SBATCH -o sor_experiment_noise_0_missing_0.3_lasso_benchmark_output.txt
#SBATCH -e sor_experiment_noise_0_missing_0.3_lasso_benchmark_error.txt
#SBATCH -J galasso_multiple_imputation
#SBATCH -p icb_cpu
#SBATCH -c 10
#SBATCH --mem=30G
#SBATCH -t 24:00:00

R CMD BATCH '--args input="/home/icb/juan.henao/lasso_benchmark/input_data/090322_single_omics_benchmark.RData" cores=10  bench_params=c(0.3,0)' single_omics_reduced_kimono_runner.R
