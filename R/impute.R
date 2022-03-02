#library(mice)
library(impute)

#setwd("/home/icb/anastasiia.grekova/")
#load("ECCB_2022/data/220221_benchmark_v2.RData")

load("data/220221_benchmark_v2.RData")

# Which data is missing?
message(paste0("# of missing values by genes ",sum(is.na(input_data$gene))))
message(paste0("# of missing values by methylation ",sum(is.na(input_data$methylation))))
message(paste0("# of missing values by cnv ",sum(is.na(input_data$cnv))))
message(paste0("# of missing values by phenotype overall ",sum(is.na(input_data$phenotype))))
message(paste0("# of missing values by phenotype: time ",sum(is.na(input_data$phenotype$srv.time))))
message(paste0("# of missing values by phenotype: status ",sum(is.na(input_data$phenotype$srv.status))))
message(paste0("# of missing values by mutation ",sum(is.na(input_data$mutation))))

THRESHOLD_MISSING <- 0.2
SEED <- 42
K <- 10

# GENES
# rows are sample IDs, columns are gene IDs
gene <- as.data.frame(input_data$gene)
# transform
gene.t <- t(gene)
print(dim(gene.t))
message("genes x samples before the knn imputation ")
geneImp <- impute.knn(gene.t, k=K, rowmax = THRESHOLD_MISSING, colmax = THRESHOLD_MISSING, rng.seed = SEED)
imputedGene <- t(geneImp$data)
message(paste0("# of missing values by genes after imputation ",sum(is.na(imputedGene))))

message(paste0("# of missing values by genes ",sum(is.na(imputed_data$gene))))
message(paste0("# of missing values by methylation ",sum(is.na(imputed_data$methylation))))
message(paste0("# of missing values by cnv ",sum(is.na(imputed_data$cnv))))
message(paste0("# of missing values by phenotype overall ",sum(is.na(imputed_data$phenotype))))
message(paste0("# of missing values by phenotype: time ",sum(is.na(imputed_data$phenotype$srv.time))))
message(paste0("# of missing values by phenotype: status ",sum(is.na(imputed_data$phenotype$srv.status))))
message(paste0("# of missing values by mutation ",sum(is.na(imputed_data$mutation))))

#geneTidy <- gene[which(rowMeans(!is.na(gene)) > 1-THRESHOLD_MISSING), which(colMeans(!is.na(gene)) > 1-THRESHOLD_MISSING)]
#print(dim(geneTidy))
#message("genes after removing missing values ")
#geneImp <- mice(geneTidy, meth="rf", ntree=10)
#geneTidyImp <- complete(geneImp)

#METHYLATION
# rows are sample IDs, columns are gene methylation sites
methylation <- as.data.frame(input_data$methylation)
# transform
methylation.t <- t(methylation)
print(dim(methylation.t))
message("methylation x samples before the knn imputation")
methylationImp <- impute.knn(methylation.t, k=K, rowmax = THRESHOLD_MISSING, colmax = THRESHOLD_MISSING, rng.seed = SEED)
imputedMethylation <- t(methylationImp$data)
message(paste0("# of missing values by methylation after imputation ",sum(is.na(imputedMethylation))))

imputed_data <- input_data
imputed_data$gene <- as.data.table(imputedGene)
imputed_data$methylation <- as.data.table(imputedMethylation)


#methylationTidy <- methylation[which(rowMeans(!is.na(methylation)) > 1-THRESHOLD_MISSING), which(colMeans(!is.na(methylation)) > 1-THRESHOLD_MISSING)]
#print(dim(methylationTidy))
#message("methylations after removing missing values ")
#methylationImp <- mice(methylationTidy, meth="rf", ntree=10)
#methylationTidyImp <- complete(methylationImp)

#save.image(file = "workspace_output.RData")
