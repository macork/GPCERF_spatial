#' @title
#' Calculate Posterior Standard Deviations for nnGP Model
#'
#' @description
#' Calculates the posterior standard deviation of a point on the CERF based on
#' the nnGP model.
#'
#' @param hyperparam The values of hyperparameters in the GP model.
#' @param w  The exposure level for the point of interest on the CERF.
#' @param GPS_w The GPS for all samples when their exposure levels are set
#' at \code{w}.
#' @param obs_ord A matrix of two columns. The first column is the observed
#' exposure levels of all samples; the second is the GPS at the observed
#' exposure levels. The rows are in ascending order for the first column.
#' @param sigma2 A scaler representing \code{sigma^2}.
#' @param n_neighbor Number of nearest neighbors on one side
#' (see also \code{expand}).
#' @param expand A scaling factor to determine the total number of nearest
#' neighbors. The total is \code{2*expand*n_neighbor}.
#'
#' @return
#' The posterior standard deviation of the estimated CERF at \code{w}.
#'
#' @keywords internal
#'
compute_posterior_sd_nn <-  function(hyperparam,
                                     w,
                                     GPS_w,
                                     obs_ord,
                                     sigma2,
                                     n_neighbor = 10,
                                     expand = 1){

  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  n <- length(GPS_w)

  if(w >= obs_ord[nrow(obs_ord),1]){
    idx_all <- seq( nrow(obs_ord) - expand*n_neighbor + 1, nrow(obs_ord), 1)
  }else{
    idx_anchor <- which.max(obs_ord[,1]>=w)
    idx_start <- max(1, idx_anchor - n_neighbor*expand)
    idx_end <- min(nrow(obs_ord), idx_anchor + n_neighbor*expand)
    if(idx_end == nrow(obs_ord)){
      idx_all <- seq(idx_end - n_neighbor*2*expand + 1, idx_end, 1)
    }else{
      idx_all <- seq(idx_start, idx_start+n_neighbor*2*expand-1, 1)
    }
  }

  obs_use <- t(t(obs_ord[idx_all,])*(1/sqrt(c(alpha, beta))))
  cov_use_inv <- chol2inv(chol(sigma2*(g_sigma*exp(-as.matrix(dist(obs_use))^2) +
                                           diag(nrow(obs_use)))))
  obs_new <- t(t(cbind(w, GPS_w))*(1/sqrt(c(alpha, beta))))

  #within variance
  sigma_sq1 <- (1+g_sigma)*sigma2/n

  #cross variance
  cross_cov <- sigma2*g_sigma*exp(-spatstat.geom::crossdist(obs_new[,1],obs_new[,2],
                                                                obs_use[,1],obs_use[,2])^2)

  sigma_sq2 <- c(calc_cross(cross_cov, cov_use_inv))/n^2
  posterior_sd <- sqrt(sigma_sq1 - sigma_sq2 + sigma2)

  logger::log_debug("w: {w}, sigma_sq1: {sigma_sq1}, sigma_sq2: {sigma_sq2},",
                    "sigma2: {sigma2}")



  return(posterior_sd)
}
