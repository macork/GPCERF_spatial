#' @title
#' Title
#'
#' @description
#'
#' @param param
#' @param sim.data
#' @param GPS
#'
#' @return
#' @export
#'
#' @examples
noise.est = function(param, sim.data, GPS){
  # browser()
  n.sample = nrow(sim.data)
  w.all = rep(sim.data$treat, 2)
  GPS.all = rep(GPS, 2)
  Sigma.all = param[3]*exp(-as.matrix(dist(cbind(w.all*sqrt(param[1]), GPS.all*sqrt(param[2]))))) +
    diag(2*n.sample)

  sd(sim.data$Y - Sigma.all[1:n.sample,-(1:n.sample)]%*%chol2inv(chol(Sigma.all[-(1:n.sample),-(1:n.sample)]))%*%sim.data$Y)
}
