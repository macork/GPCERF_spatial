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
  w.mean = sum(w*weights)
  w.sd = sqrt(sum((w-w.mean)^2*weights))
  w.trans = (w-w.mean)/w.sd

  X.mean = colSums(X*weights)
  X.cov = (t(X) - X.mean)%*%diag(weights)%*%t(t(X)-X.mean)
  X.trans = t(t(solve(chol(X.cov)))%*%(t(X)-X.mean))

  c(w.trans%*%diag(weights)%*%X.trans)
}
