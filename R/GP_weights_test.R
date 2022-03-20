#' @title
#' Calculate Weights for Estimation of a Point on CERF
#'
#' @description
#' Calculates the weights of observed outcomes which is then used to estimate
#' the posterior mean of CERF at a given exposure level.
#'
#' @param w param's A scalar of exposure level of interest.
#' @param w.obs A vector of observed exposure levels of all samples.
#' @param obs.use A matrix of two columns. First column is the observed exposure levels of all
#' samples; second is the GPS at the observed exposure levels.
#' @param param A vector of hyperparameters for the GP. (alpha, beta, gamma)
#' @param inv.Sigma.obs Inverse of the covariance matrix between observed samples.
#' @param e_gps_pred A vector of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param kernel.fn The covariance function of GP.
#'
#' @return
#' A vector of the weights assigned to each sample for the calculate of posterior mean
#' of CERF at \code{w}.
#' @export
#'
#' @examples
GP.weights.test = function(w, w.obs, obs.use, param, inv.Sigma.obs,
                           e_gps_pred, e_gps_std,
                           kernel.fn = function(x) exp(-x^2)){

  # param[1]: alpha, param[2]: beta, param[3]: gamma
  GPS.new = stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)

  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1],
                                                            obs.new[,2],
                                                            obs.use[,1],
                                                            obs.use[,2]))
  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  # system.time(eigenMapMatMult(Sigma.cross, inv.Sigma.obs))
  c((rep(1/length(w.obs),length(w.obs))%*%Sigma.cross)%*%inv.Sigma.obs)
}
