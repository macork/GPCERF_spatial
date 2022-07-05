#' @title
#' Calculate Covariate Balance
#'
#' @description
#' Calculates the weighted correlation between a list of covariates and an
#' exposure. The covariance function of the GP defines weights.
#'
#' @param w A vector of exposure values across all subjects.
#' @param X A matrix of covariate values. Subjects in rows and covariates in columns.
#' @param weights A vector of weights assigned to all subjects based on the trained GP.
#'
#' @return
#' A vector of correlations between w and each column of X.
#'
#' @keywords internal
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
