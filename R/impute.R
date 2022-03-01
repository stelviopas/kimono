library(mice)
setwd("/home/anastasiia.grekova/kimono/")
load("data/220221_benchmark.RData")

# Which data is missing?
message(paste0("# of missing values by genes ",sum(is.na(input_data$gene))))
message(paste0("# of missing values by methylation ",sum(is.na(input_data$methylation))))
message(paste0("# of missing values by cnv ",sum(is.na(input_data$cnv))))
message(paste0("# of missing values by phenotype overall ",sum(is.na(input_data$phenotype))))
message(paste0("# of missing values by phenotype: time ",sum(is.na(input_data$phenotype$srv.time))))
message(paste0("# of missing values by phenotype: status ",sum(is.na(input_data$phenotype$srv.status))))
message(paste0("# of missing values by mutation ",sum(is.na(input_data$mutation))))

THRESHOLD_MISSING <- 0.2

# GENES
# rows are sample IDs, columns are gene IDs
gene <- as.data.frame(input_data$gene)
print(dim(gene))
message("genes before removing missing values ")
geneTidy <- gene[which(rowMeans(!is.na(gene)) > 1-THRESHOLD_MISSING), which(colMeans(!is.na(gene)) > 1-THRESHOLD_MISSING)]
print(dim(geneTidy))
message("genes after removing missing values ")
geneImp <- mice(geneTidy, meth="rf", ntree=10)
geneTidyImp <- complete(geneImp)
# METHYLATION
# rows are sample IDs, columns are gene methylation sites
methylation <- as.data.frame(input_data$methylation)
print(dim(methylation))
message("methylations before removing missing values ")
methylationTidy <- methylation[which(rowMeans(!is.na(methylation)) > 1-THRESHOLD_MISSING), which(colMeans(!is.na(methylation)) > 1-THRESHOLD_MISSING)]
print(dim(methylationTidy))
message("methylations after removing missing values ")
methylationImp <- mice(methylationTidy, meth="rf", ntree=10)
methylationTidyImp <- complete(methylationImp)
save.image(file = "workspace_output.RData")
