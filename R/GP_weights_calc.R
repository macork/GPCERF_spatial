#' @title
#' Title
#'
#' @param w
#' @param w.obs
#' @param obs.use
#' @param param
#' @param inv.Sigma.obs
#' @param e_gps_pred
#' @param e_gps_std
#' @param kernel.fn
#'
#' @return
#' @export
#'
#' @examples
GP.weights.calc = function(w, w.obs, obs.use, param, inv.Sigma.obs, e_gps_pred, e_gps_std,
                           kernel.fn = function(x) exp(-x^2)){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std)

  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                                            obs.use[,1], obs.use[,2]))
  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  weights.all = Sigma.cross%*%inv.Sigma.obs
  weights.all
}
