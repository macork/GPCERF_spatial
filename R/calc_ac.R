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
calc.ac = function(w, X, weights){
  w.mean = sum(w*weights)
  w.sd = sqrt(sum((w-w.mean)^2*weights))
  w.trans = (w-w.mean)/w.sd

  X.mean = colSums(X*weights)
  X.cov = (t(X) - X.mean)%*%diag(weights)%*%t(t(X)-X.mean)
  X.trans = t(t(solve(chol(X.cov)))%*%(t(X)-X.mean))

  c(w.trans%*%diag(weights)%*%X.trans)
}
