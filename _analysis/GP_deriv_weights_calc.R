#' @title
#' Calculate Derivatives of CERF
#'
#' @description
#' Calculate the weights assigned to each observed outcome when deriving the
#' posterior mean of the first derivative of CERF at a given exposure level.
#'
#' @param w A scalar of exposure level of interest.
#' @param w.obs A vector of observed exposure levels of all samples.
#' @param GPS.obs A vector of GPS for all samples at the observed levels of exposure.
#' @param param A vector of hyper-parameters in the GP model.
#' @param e_gps_pred A vecotor of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param kernel.fn The covariance function.
#' @param kernel.deriv.fn The partial derivative of the covariance function.
#'
#' @return
#' A vector of weights for all samples, based on which the posterior mean of the derivative of CERF at the
#' exposure level of interest is calculated.
#' @export
#'
#' @examples
GP.deriv.weights.calc = function(w, w.obs, GPS.obs, param, e_gps_pred, e_gps_std,
                                 kernel.fn = function(x) exp(-x),
                                 kernel.deriv.fn = function(x) -exp(-x)){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)
  n = length(GPS.new)

  obs.use = cbind( w.obs*sqrt(param[1]), GPS.obs*sqrt(param[2]) )
  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))^2) + diag(nrow(obs.use))
  cross.dist = spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                        obs.use[,1], obs.use[,2])
  Sigma.cross = param[3]*sqrt(param[1])*kernel.deriv.fn(cross.dist^2)*
    (2*outer(rep(w,n), w.obs, "-"))
  weights.all = Sigma.cross%*%chol2inv(chol(Sigma.obs))
  # weights.all[weights.all<0] = 0
  # weights = colMeans(weights.all)
  # weights/sum(weights)
  colMeans(weights.all)
}
