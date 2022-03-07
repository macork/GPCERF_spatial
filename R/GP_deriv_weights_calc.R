# weights for derivatives
# we switch back to dnorm to calc GPS
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
                                 kernel.fn = function(x) exp(-x^2),
                                 kernel.deriv.fn = function(x, mu, sigma) (x-mu)/sigma^2){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std)

  obs.use = cbind( w.obs*sqrt(param[1]), GPS.obs*sqrt(param[2]) )
  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use))
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                                            obs.use[,1], obs.use[,2]))*
    param[1]*outer(w, w.obs, "-") - param[2]*outer(GPS.new, GPS.obs, "-")*
    GPS.new*kernel.deriv.fn(w, e_gps_pred, e_gps_std)
  weights.all = Sigma.cross%*%chol2inv(chol(Sigma.obs))
  weights.all
}
