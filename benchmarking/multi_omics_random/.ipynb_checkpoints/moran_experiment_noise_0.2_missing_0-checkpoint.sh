#!/bin/bash

#SBATCH -o moran_experiment_noise_0.2_missing_0_lasso_benchmark_output.txt
#SBATCH -e moran_experiment_noise_0.2_missing_0_lasso_benchmark_error.txt
#SBATCH -J galasso_multiple_imputation
#SBATCH -p icb_cpu
#SBATCH -c 30
#SBATCH --mem=30G
#SBATCH -t 24:00:00

R CMD BATCH '--args input="/home/icb/juan.henao/lasso_benchmark/multi_omics_random/240222_multi_omics_benchmark.RData" cores=30  bench_params=c(0,0.2)' multi_omics_random_kimono_runner.R
