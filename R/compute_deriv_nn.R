#' @title
#' Calculate Derivatives of CERF for nnGP
#'
#' @description
#' Calculates the posterior mean of the derivative of CERF at a given
#' exposure level with nnGP.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param y_obs A vector of observed outcome values.
#' @param hyperparam A vector of hyper-parameters in the GP model.
#' @param n_neighbor Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n_neighbor}.
#' @param block_size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block_size} is faster, but requires more memory.
#' @param kernel_fn The covariance function. The input is the square of Euclidean distance.
#' @param kernel_deriv_fn The partial derivative of the covariance function. The input is the square of Euclidean distance.
#'
#' @return
#' A scalar of estimated derivative of CERF at \code{w} in nnGP.
#' @export
#'
#' @examples
#'
#' set.seed(365)
#' data <- generate_synthetic_data(sample_size = 200)
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' wi <- 4.8
#'
#' deriv_val <- compute_deriv_nn(w = wi,
#'                               w_obs = data$treat,
#'                               GPS_m = GPS_m,
#'                               y_obs = data$Y,
#'                               hyperparam = c(0.1,0.2,1),
#'                               n_neighbor = 20,
#'                               expand = 1,
#'                               block_size = 1000)
#'
compute_deriv_nn <- function(w,
                             w_obs,
                             GPS_m,
                             y_obs,
                             hyperparam,
                             n_neighbor,
                             expand,
                             block_size,
                             kernel_fn = function(x) exp(-x),
                             kernel_deriv_fn = function(x) -exp(-x)){


  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  GPS <- GPS_m$GPS
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std


  # params[1]: alpha, params[2]: beta, params[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS_w = dnorm(w, mean = e_gps_pred, sd = e_gps_std)

  n = length(GPS_w)
  n.block = ceiling(n/block_size)
  obs.raw = cbind(w_obs, GPS)
  obs.ord = obs.raw[order(obs.raw[,1]),]
  y_obs.ord = y_obs[order(obs.raw[,1])]
  #params: length 3, first scale for w, second scale for GPS,
  #third scale for exp fn
  if(w >= obs.ord[nrow(obs.ord),1]){
    idx.all = seq( nrow(obs.ord) - expand*n_neighbor + 1, nrow(obs.ord), 1)
  }else{
    idx.anchor = which.max(obs.ord[,1]>=w)
    idx.start = max(1, idx.anchor - n_neighbor*expand)
    idx.end = min(nrow(obs.ord), idx.anchor + n_neighbor*expand)
    if(idx.end == nrow(obs.ord)){
      idx.all = seq(idx.end - n_neighbor*2*expand + 1, idx.end, 1)
    }else{
      idx.all = seq(idx.start, idx.start+n_neighbor*2*expand-1, 1)
    }
  }

  obs.use = t(t(obs.ord[idx.all,])*(1/sqrt(c(alpha, beta))))
  y.use = y_obs.ord[idx.all]

  obs.new = t(t(cbind(w, GPS_w))*(1/sqrt(c(alpha, beta))))
  id.all = split(1:n, ceiling(seq_along(1:n)/n.block))
  Sigma.obs = g_sigma*kernel_fn(as.matrix(dist(obs.use))^2) + diag(nrow(obs.use))
  Sigma.obs.inv = chol2inv(chol(Sigma.obs))

  all.weights = sapply(id.all, function(id.ind){
    cross.dist = spatstat.geom::crossdist(obs.new[id.ind,1], obs.new[id.ind,2],
                                          obs.use[,1], obs.use[,2])
    Sigma.cross = g_sigma*(1/alpha)*(-2*outer(rep(w,length(id.ind))*(1/alpha), obs.use[,1], "-"))*
      kernel_deriv_fn(cross.dist^2)
    #mean
    wght = Sigma.cross%*%Sigma.obs.inv
    wght[wght<0] = 0
    colSums(wght)
  })
  weights = rowSums(all.weights)/n
  weights = weights/sum(weights)
  weights%*%y.use
}
