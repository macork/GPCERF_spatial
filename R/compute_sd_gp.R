#' @title
#' Compute Posterior Credible Interval
#'
#' @description
#' Computes posterior credible interval for requested exposure level.
#'
#' @param w A scalar of exposure level of interest.
#' @param scaled_obs A matrix of two columns.
#'   - First column is the scaled GPS value of all samples (GPS * 1/sqrt(alpha))
#'   - Second column is the scaled exposure value of all samples (w * 1/sqrt(beta))
#' @param hyperparam A vector of hyper-parameters for the GP.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: gamma/sigma
#' @param sigma  A scaler that represents noise.
#' @param GPS_m  A data.table of GPS vectors.
#'   - Column 1: A vector of estimated GPS evaluated at the observed exposure levels.
#'   - Column 2: Estimated conditional means of the exposure given covariates
#'               for all samples (e_gps_pred).
#'   - Column 3: Estimated conditional standard deviation of the exposure given
#'               covariates for all samples (e_gps_std).
#' @param kernel_fn The covariance function of GP.
#'
#' @keywords internal
#'
#' @return
#' Posterior credible interval (scaler) for the requested exposure level (w).
#'

compute_sd_gp <- function(w,
                          scaled_obs,
                          hyperparam,
                          sigma,
                          GPS_m,
                          kernel_fn = function(x) exp(-x^2)){


  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]

  n <- nrow(scaled_obs)

  # Compute GPS for requested w
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std
  GPS_w <- stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)

  # Compute helper matrix for the new w and corresponding GPS.
  scaled_w <- cbind( w/sqrt(1/alpha), GPS_w/sqrt(1/beta))

  scaled_combined <- rbind(scaled_w, scaled_obs)
  Sigma_all <- (g_sigma*kernel_fn(as.matrix(stats::dist(scaled_combined))) + diag(n*2))*sigma^2
  Sigma_within_w <- Sigma_all[1:n, 1:n]
  Sigma_cross <- Sigma_all[1:n, -(1:n)]
  Sigma_within_obs <- Sigma_all[-(1:n), -(1:n)]
  Sigma_conditional <- Sigma_within_w - Sigma_cross%*%chol2inv(chol(Sigma_within_obs))%*%t(Sigma_cross)
  posterior_sd <- sqrt(rep(1/n,n)%*%Sigma_conditional%*%rep(1/n,n) + sigma^2)
  return(posterior_sd)
}
