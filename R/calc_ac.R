#' @title
#' Calculate Covariate Balance
#'
#' @description
#' Calculate weighted correlation between a list of covariates and an exposure.
#' Weights are defined by the covariance function of the GP.
#'
#' @param w A vector exposure values across all subjects.
#' @param X A matrix of covariate values. Subjects in rows and covariates in columns.
#' @param weights A vector of weights assigned to all subjects based on the trained GP.
#'
#' @return
#' A vector of correlations between w and each column of X.
#' @export
#'
#' @examples
#'
#' set.seed(429)
#'
#' # generate data
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#'
#' # generate random weights
#' weights <- runif(nrow(data))
#' weights <- weights/sum(weights)
#'
#' # covariate matrix
#' design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
#'
#' cb <- calc_ac(w = data$treat, X = design_mt, weights=weights)
#'
calc_ac <- function(w, X, weights){

  w_mean <- sum(w*weights)
  w_sd <- sqrt(sum((w-w_mean)^2*weights))
  w_trans <- (w-w_mean)/w_sd

  X_mean <- colSums(X*weights)
  X_cov <- (t(X) - X_mean)%*%diag(weights)%*%t(t(X)-X_mean)
  X_trans <- t(t(solve(chol(X_cov)))%*%(t(X)-X_mean))

  return(c(w_trans%*%diag(weights)%*%X_trans))
}
