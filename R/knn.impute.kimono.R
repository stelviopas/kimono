knn.impute.kimono <- function(input_data, prior_network, seed, core = 2, infer_missing_prior = TRUE){

  library(impute)
  library(data.table)
  library(kimono)

  message("KNN imputation + KiMONo")
  imputed_data <- input_data
  # Which data is missing?
  # "gene", "methylation", "cnv", "phenotype"
  for (data.table.name in names(input_data)) {
    mis.val.count <- sum(is.na(input_data[[data.table.name]]))
    message(paste0("# of missing values by ", data.table.name, ": ", mis.val.count))
  }
  message("Starting KNN imputation ...")
  # impute missing values on each level (e.g. gene or cnv data.table) separately
  for (data.table.name in names(input_data)) {
    data.table <- input_data[[data.table.name]]
    mis.val.count <- sum(is.na(data.table))
    if (mis.val.count == 0) {
      message(paste0("No missing data. Skip ", data.table.name))
    }
    else{
      message(paste0("Imputing ", data.table.name))
      data <- as.data.frame(data.table)
      # transform
      # originally, knn.impute requires rows columns as samples
      # this is important since default behavior is any row(e.g. gene)
      # with more than rowmax% (default: 50%) missing are imputed using the overall mean per sample (column)
      # additionally, if any column has more than colmax% (default: 80%) missing data, the program halts and reports an error
      data.t <- t(data)
      dataImp <- impute.knn(data.t, rng.seed = seed)
      dataImp <- t(dataImp$data)
      message("Finished")
      imputed_data[[data.table.name]] <- as.data.table(dataImp)
    }
  }
  message("Calling KiMONo ...")
  return(kimono(imputed_data, prior_network,
                core = core,
                infer_missing_prior = infer_missing_prior))
}
