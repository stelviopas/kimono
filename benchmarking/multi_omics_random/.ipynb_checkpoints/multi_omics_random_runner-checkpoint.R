var <- paste0("#!/bin/bash

#SBATCH -o moran_experiment_noise_jpivot_missing_ipivot_lasso_benchmark_output.txt
#SBATCH -e moran_experiment_noise_jpivot_missing_ipivot_lasso_benchmark_error.txt
#SBATCH -J galasso_multiple_imputation
#SBATCH -p icb_cpu
#SBATCH -c 10
#SBATCH --mem=30G
#SBATCH -t 24:00:00

R CMD BATCH '--args input=","\"/home/icb/juan.henao/lasso_benchmark/input_data/090322_multi_omics_benchmark.RData\""," cores=10  bench_params=white' multi_omics_random_kimono_runner.R")

# c(i,j)

perc_miss <- c(0.0,0.10,0.30,0.50,0.70,0.90)
noise_per <- c(0.0,1.0,1.3,1.5,1.7,2.0)

#perc_miss <- c(0.1)
#noise_per <- c(0)

for(i in perc_miss){
for(j in noise_per){
    repl  <- paste0("c(",i,",",j,")")
tmp = gsub(pattern = "white",replacement = repl,x = var)
tmp <- gsub(pattern = "jpivot",replacement = j,x = tmp)
tmp <- gsub(pattern = "ipivot",replacement = i,x = tmp)
    #print(tmp)
file_name <- paste0("moran_experiment_noise_",j,"_missing_",i,".sh")
    file_con <- file(file_name)
writeLines(text = c(tmp),file_con)
    close(file_con)
    #write.table(tmp,file_name,quote = F,row.names = F,col.names = F)
system(command=paste0("bash ./",file_name))
}
}