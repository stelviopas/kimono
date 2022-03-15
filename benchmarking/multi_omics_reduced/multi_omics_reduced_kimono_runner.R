library(kimono)

args = (commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    input = "/home/icb/juan.henao/lasso_benchmark/input_data/090322_multi_omics_benchmark.RData"
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
    
    if(missingness != 0){
        gender.two <- which(input_data[['phenotype']]$gender == 2)
    to.deletion <- sample(x = nrow(input_data[['phenotype']]), size = round(nrow(input_data[['phenotype']]) * missingness))
    to.deletion <- ifelse(to.deletion %in% gender.two,yes = to.deletion + 1,no = to.deletion)
    all.pos <- which(!rownames(input_data[['phenotype']]) %in% to.deletion)  

    input_data[['gene']] <- input_data[['gene']][all.pos,]
    input_data[['cnv']] <- input_data[['cnv']][all.pos,]
    input_data[['methylation']] <- input_data[['methylation']][all.pos,]
    input_data[['phenotype']] <- input_data[['phenotype']][all.pos,]
    
    for(layer in c('gene','cnv','methylation')){
    layer_scaled <- scale(input_data[[layer]])
        col_removes <- apply(layer_scaled , 2, stats::var) == 0
        input_data[[layer]] <- input_data[[layer]][,!col_removes, with=FALSE]
    }    
    }
    
    if(white_noise != 0){
    for(layer in c("gene","cnv","methylation")){
        sds_matrix <- matrix(0,nrow = nrow(input_data[[layer]]),ncol=ncol(input_data[[layer]]))
    sds <- apply(as.matrix(input_data[[layer]]),2,sd)
    for(x in seq(length(sds))){
    sds_matrix[,x] <- rnorm(n = nrow(input_data[[layer]]),mean = 0, sd = sds[x] * white_noise)
    }
    input_data[[layer]] <- input_data[[layer]] + sds_matrix
    
    input_data[[layer]] <- setDT(input_data[[layer]])
    }
    pre.methylation <- as.matrix(input_data[['methylation']])
    pre.methylation[which(pre.methylation < 0)] <- 0
    pre.methylation[which(pre.methylation > 1)] <- 1
    input_data[['methylation']] <- setDT(as.data.frame(pre.methylation))
    }
    
    
    
    
    network <- kimono(input_data, prior_network ,core = cores, infer_missing_prior = FALSE)
    write.csv(network,paste0("mor_experiment_noise_",white_noise,"_missing_",missingness,"_iteration_",iter,".csv"), quote = FALSE, row.names = FALSE)
}

