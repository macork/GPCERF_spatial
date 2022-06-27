#' @title
#' Calculate Posterior Standard Deviations for nnGP Model
#'
#' @description
#' Calculates the posterior standard deviation of a point on the CERF based on the nnGP model.
#'
#' @param hyperparam Values of hyperparameters in the GP model.
#' @param w  The exposure level for the point of interest on the CERF.
#' @param GPS_w The GPS for all samples when their exposure levels are set at \code{w}.
#' @param obs_ord A matrix of two columns. First column is the observed exposure levels of all
#' samples; second is the GPS at the observed exposure levels. The rows are in ascending order
#' for the first column.
#' @param sigma2 A scaler representing \code{sigma^2}
#' @param n_neighbor Number of nearest neighbors on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbors. The total is \code{2*expand*n_neighbor}.
#'
#' @return
#' The posterior standard deviation of the estimated CERF at \code{w}.
#' @export
#'
#' @examples
#'
#' set.seed(3089)
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' # Hyperparameter
#' hyperparam <- c(0.1, 0.2, 1)
#' n_neighbor <- 10
#' expand <- 1
#' block_size <- 10000
#'
#' # Exposure level
#' wi <- 0.4
#'
#' # Estimate GPS for the exposure level
#' GPS_w = dnorm(wi,
#'               mean = GPS_m$e_gps_pred,
#'               sd = GPS_m$e_gps_std, log = TRUE)
#'
#' # Order data for easy selection
#' coord_obs = cbind(data$treat, GPS_m$GPS)
#' y_use <- data$Y
#'
#' obs_ord <- coord_obs[order(coord_obs[,1]),]
#' y_use_ord <- y_use[order(coord_obs[,1])]
#'
#' # compute noise
#' noise <- estimate_noise_nn(hyperparam = hyperparam,
#'                            w_obs = data$treat,
#'                            GPS_obs = GPS_m$GPS,
#'                            y_obs = y_use_ord,
#'                            n_neighbor = n_neighbor)
#'
#' # compute posterior standard deviation
#' pst_sd <- compute_posterior_sd_nn(hyperparam = hyperparam,
#'                                   w = wi,
#'                                   GPS_w = GPS_w,
#'                                   obs_ord = obs_ord,
#'                                   sigma2 = noise,
#'                                   n_neighbor = 20,
#'                                   expand = 1)
#'
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
