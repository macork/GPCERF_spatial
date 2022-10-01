#' @title
#' Compute Weighted Correlation
#'
#' @description
#' Computes weighted correlation of the observational data based on weights
#' achieved by Gaussian Process.
#'
#' @param data A data.table of observational data with the following colums:
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#' @param weights A vector of weights for each observation data.
#'
#' @return
#' A vector of covariate balance
#'
#' @export
#'
#' @examples
#'
#' set.seed(124)
#' mydata <- generate_synthetic_data(sample_size = 200)
#' data.table::setDT(mydata)
#' weights <- runif(nrow(mydata))
#' compute_w_corr(mydata, weights)
#'
compute_w_corr <- function(data, weights){


  if (!is.data.table(data)){
    stop(paste0("The data should be a data.table. ",
                "Current format: ", class(data)[1]))
  }

  if (nrow(data) != length(weights)){
    stop(paste0("Number of data samples (", nrow(data), ") and length of ",
                "weights (", length(weights),") should be equal."))
  }

  # TODO: model.matrix will create dummy vairables for factors.
  # Double check with Boyu.
  conf_names <- colnames(data[,3:ncol(data)])
  frml <- paste("~",paste(conf_names, collapse = "+"), "-1", sep = "")

  w_obs <- data[[2]]

  x_design <- model.matrix(as.formula(frml), data = data)
  w_mean <- sum(w_obs*weights)
  w_sd <- sqrt(sum((w_obs - w_mean)^2*weights))
  w_stan <- (w_obs - w_mean)/w_sd

  x_mean <- colSums(x_design*weights)
  x_cov <- (t(x_design) - x_mean)%*%diag(weights)%*%t(t(x_design) - x_mean)
  x_cov_ev <- eigen(x_cov)$values
  # when x_cov is rank deficient, return NA for all covariate balance
  if(sum(x_cov_ev<=0)==0){
    x_stan <- t(t(solve(chol(x_cov)))%*%(t(x_design) - x_mean))
    covariate_balance <- abs(c(t(x_stan)%*%diag(weights)%*%w_stan))
  }else{
    covariate_balance = rep(NA, nrow(x_cov))
  }
  return(covariate_balance)
}
