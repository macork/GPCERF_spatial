# weights for derivatives
# we switch back to dnorm to calc GPS
#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param w param's description
#' @param w.obs param's description
#' @param GPS.obs param's description
#' @param param param's description
#' @param e_gps_pred param's description
#' @param e_gps_std param's description
#' @param kernel.fn param's description
#' @param kernel.deriv.fn param's description
#'
#' @return
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
