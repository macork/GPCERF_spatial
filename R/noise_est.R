#' @title
#' Estimate the Standard Deviation of the Nugget Term in Full Gaussian Process
#'
#' @description
#' Estimates the standard deviations of the nugget term in full GP by calculating
#' the standard deviations of the residuals.
#'
#' @param param A vector of hyper-parameter values for the full GP.
#' @param data A data frame containing outcome, exposure and all covariates.
#' @param GPS A vector of estimated GPS at the observed exposure levels.
#'
#' @return
#' A scalar of estimated standard deviation of the nugget term in full GP.
#' @export
#'
#' @examples
#'
noise.est = function(param, data, GPS){

  n.sample = nrow(sim.data)
  w.all = rep(data$treat, 2)
  GPS.all = rep(GPS, 2)
  Sigma.all = param[3]*exp(-as.matrix(dist(cbind(w.all*sqrt(param[1]),
                                                 GPS.all*sqrt(param[2]))))) +
    diag(2*n.sample)

  sd(data$Y - Sigma.all[1:n.sample,
                        -(1:n.sample)]%*%chol2inv(
                          chol(Sigma.all[-(1:n.sample),
                                         -(1:n.sample)]))%*%sim.data$Y)
}
