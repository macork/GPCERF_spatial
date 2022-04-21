#' @title
#' Calculate Derivatives of CERF for nnGP
#'
#' @description
#' Calculate the posterior mean of the derivative of CERF at a given
#' exposure level with nnGP.
#'
#' @param w A scalar of exposure level of interest.
#' @param w.obs A vector of observed exposure levels of all samples.
#' @param GPS.obs A vector of GPS for all samples at the observed levels of exposure.
#' @param y.obs A vector of observed outcome values.
#' @param params A vector of hyper-parameters in the GP model.
#' @param e_gps_pred A vecotor of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param n.neighbour Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n.neighbour}.
#' @param block.size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block.size} is faster, but requires more memory.
#' @param kernel.fn The covariance function. The input is the square of Euclidean distance.
#' @param kernel.deriv.fn The partial derivative of the covariance function. The input is the square of Euclidean distance.
#'
#' @return
#' A scalar of estimated derivative of CERF at \code{w} in nnGP.
#' @export
#'
#' @examples
deriv.nn.fast = function(w, w.obs, GPS.obs, y.obs, params, e_gps_pred, e_gps_std,
                                 n.neighbour, expand, block.size,
                                 kernel.fn = function(x) exp(-x),
                                 kernel.deriv.fn = function(x) -exp(-x)){
  # params[1]: alpha, params[2]: beta, params[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)

  n = length(GPS.new)
  n.block = ceiling(n/block.size)
  obs.raw = cbind(w.obs, GPS.obs)
  obs.ord = obs.raw[order(obs.raw[,1]),]
  y.obs.ord = y.obs[order(obs.raw[,1])]
  #params: length 3, first scale for w, second scale for GPS,
  #third scale for exp fn
  if(w >= obs.ord[nrow(obs.ord),1]){
    idx.all = seq( nrow(obs.ord) - expand*n.neighbour + 1, nrow(obs.ord), 1)
  }else{
    idx.anchor = which.max(obs.ord[,1]>=w)
    idx.start = max(1, idx.anchor - n.neighbour*expand)
    idx.end = min(nrow(obs.ord), idx.anchor + n.neighbour*expand)
    if(idx.end == nrow(obs.ord)){
      idx.all = seq(idx.end - n.neighbour*2*expand + 1, idx.end, 1)
    }else{
      idx.all = seq(idx.start, idx.start+n.neighbour*2*expand-1, 1)
    }
  }

  obs.use = t(t(obs.ord[idx.all,])*sqrt(params[1:2]))
  y.use = y.obs.ord[idx.all]

  obs.new = t(t(cbind(w, GPS.new))*sqrt(params[1:2]))
  id.all = split(1:n, ceiling(seq_along(1:n)/n.block))
  Sigma.obs = params[3]*kernel.fn(as.matrix(dist(obs.use))^2) + diag(nrow(obs.use))
  Sigma.obs.inv = chol2inv(chol(Sigma.obs))

  all.weights = sapply(id.all, function(id.ind){
    cross.dist = spatstat.geom::crossdist(obs.new[id.ind,1], obs.new[id.ind,2],
                                          obs.use[,1], obs.use[,2])
    Sigma.cross = params[3]*sqrt(params[1])*(2*outer(rep(w,length(id.ind))*params[1], obs.use[,1], "-"))*
      kernel.deriv.fn(cross.dist^2)
    #mean
    w = Sigma.cross%*%Sigma.obs.inv
    colSums(w)
  })
  weights = rowSums(all.weights)/n
  weights%*%y.use
}
