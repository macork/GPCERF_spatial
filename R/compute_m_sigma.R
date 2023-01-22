#' @title
#' Compute mean, credible interval, and covariate balance in full Gaussian
#'  process (GP)
#'
#' @description
#' Calculates the induced covariate balance associated with one hyper-parameter
#' configuration in full GP.
#'
#' @param hyperparam A vector of values of hyper-parameters.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: g_sigma (gamma / sigma)
#' @param data A  data.frame containing all data including outcome, exposure
#' and covariates. In the following order:
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#' @param w A vector of exposure levels at which the CERF is estimated.
#' @param GPS_m A data.frame of GPS vectors.
#'   - Column 1: A vector of estimated GPS evaluated at the observed exposure
#'   levels.
#'   - Column 2: Estimated conditional means of the exposure given covariates
#'               for all samples (e_gps_pred).
#'   - Column 3: Estimated conditional standard deviation of the exposure given
#'               covariates for all samples (e_gps_std).
#' @param tuning The function is used for parameter tuning (default = TRUE)
#' or estimation (FALSE)
#' @param kernel_fn The covariance function of GP.
#'
#' @return
#' A list containing two elements:
#'   - A vector of absolute weighted correlation of each covariate to the
#'   exposure, which is the metric for covariate balance
#'   - An estimated CERF at \code{w_all} based on the hyper-parameter values in
#'   \code{param}.
#'
#' @keywords internal
#'
compute_m_sigma <- function(hyperparam, data, w, GPS_m, tuning,
                            kernel_fn = function(x) exp(-x ^ 2)){

  param <- unlist(hyperparam)

  GPS <- GPS_m$GPS
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std

  # mi(w)
  # param 1: alpha
  # param 2: beta
  # param 3: ratio gamma/sigma

  alpha <- param[1]
  beta  <- param[2]
  g_sigma <- param[3]

  w_obs <- data[[2]]

  scaled_obs <- cbind(w_obs * sqrt(1 / alpha), GPS * sqrt(1 / beta))
  sigma_obs <- g_sigma * kernel_fn(as.matrix(dist(scaled_obs))) +
               diag(nrow(scaled_obs))

  inv_sigma_obs <- compute_inverse(sigma_obs)

  # Estimate noise
  if(!tuning) {
    noise_est <- estimate_noise_gp(data = data,
                                   sigma_obs, inv_sigma_obs)
  }

  col_all_list <- lapply(w,
                         function(w_instance) {

    # compute weights
    weights_res <- compute_weight_gp(w = w_instance,
                                     w_obs = w_obs,
                                     scaled_obs = scaled_obs,
                                     hyperparam = hyperparam,
                                     inv_sigma_obs = inv_sigma_obs,
                                     GPS_m = GPS_m,
                                     est_sd = !tuning,
                                     kernel_fn = kernel_fn)

    weights_final <- weights_res$weight
    weights_final[weights_final < 0] <- 0
    if(sum(weights_final) > 0) {
      weights_final <- weights_final / sum(weights_final)
    }

    # weigts.final = invers of paranthesis * kappa
    # est is the same as m in the paper.

    if(!tuning) {
      est <- data$Y %*% weights_final
      pst_sd <- noise_est * sqrt(weights_res$sd_scaled ^ 2 + 1)
    }
    else{
      est <- NA
      pst_sd <- NA
    }
    covariate_balance <- compute_w_corr(w = data[[2]],
                                        confounders = data[, 3:ncol(data)],
                                        weights_final)
    c(covariate_balance, est, pst_sd)
  })

  col_all <- do.call(cbind, col_all_list)

  n_confounders <- nrow(col_all) - 2 # est, pst_sd
  est_index <- nrow(col_all) -1
  pst_index <- nrow(col_all)

  list(cb = rowMeans(col_all[1:n_confounders, ], na.rm = T),
       est = col_all[est_index, ],
       pst = col_all[pst_index, ])
}
