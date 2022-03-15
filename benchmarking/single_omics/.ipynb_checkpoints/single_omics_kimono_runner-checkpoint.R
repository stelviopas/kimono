library(kimono)

args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    input = "/home/icb/juan.henao/lasso_benchmark/input_data/090322_single_omics_benchmark.RData"
    cores = 2
    bench_params = c(0.10,1)
}else{
    cat("using CMD args")
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

white_noise <- bench_params[2]
missingness <- bench_params[1]

for(iter in seq(5)){
    load(input)
    set.seed(iter)
    
    if(white_noise != 0){
    sds_matrix <- matrix(0,nrow = nrow(input_data[["gene"]]),ncol=ncol(input_data[["gene"]]))
    sds <- apply(as.matrix(input_data[["gene"]]),2,sd)
    for(x in seq(length(sds))){
    sds_matrix[,x] <- rnorm(n = nrow(input_data[['gene']]),mean = 0, sd = sds[x] * white_noise)
    }
    input_data[["gene"]] <- input_data[["gene"]] + sds_matrix
    
    input_data[["gene"]] <- setDT(input_data[["gene"]])
    }
    
    if(missingness != 0){
        dat <- as.matrix(input_data[["gene"]])
    n <- round(length(dat) * missingness)
    positions_sample <- sample(x = length(dat), size = n)
    dat[positions_sample] <- NA
    input_data[["gene"]] <- setDT(as.data.frame(dat))
    }
    
    
    network <- kimono(input_data, prior_network ,core = cores, infer_missing_prior = FALSE)
    write.csv(network,paste0("so_experiment_noise_",white_noise,"_missing_",missingness,"_iteration_",iter,".csv"), quote = FALSE, row.names = FALSE)
}
