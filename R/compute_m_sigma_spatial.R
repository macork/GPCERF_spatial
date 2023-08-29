#' @title
#' Compute mean, credible interval, and covariate balance in standard Gaussian
#'  process (GP)
#'
#' @description
#' Calculates the induced covariate balance associated with one hyper-parameter
#' configuration in standard GP.
#'
#' @param hyperparam A vector of values of hyper-parameters.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: g_sigma (gamma / sigma)
#' @param outcome_data A  vector of outcome data.
#' @param treatment_data A vector of treatment data.
#' @param covariates_data A data frame of covariates data.
#' @param w A vector of exposure levels at which the CERF is estimated.
#' @param gps_m An S3 gps object including:
#'   gps: A data.frame of GPS vectors.
#'     - Column 1: GPS
#'     - Column 2: Prediction of exposure for covariate of each data sample
#'     (e_gps_pred).
#'     - Column 3: Standard deviation of  e_gps (e_gps_std)
#'   used_params:
#'     - dnorm_log: TRUE or FLASE
#' @param tuning The function is used for parameter tuning (default = TRUE)
#' or estimation (FALSE)
#' @param spatial_coords Spatial coordinates at location of points used to calculate relative distance
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
compute_m_sigma_spatial <- function(hyperparam, outcome_data, treatment_data,
                            covariates_data, w, gps_m, tuning, spatial_coords,
                            kernel_fn = function(x) exp(-x ^ 2)) {

  param <- unlist(hyperparam)

  gps <- gps_m$gps$gps

  # mi(w)
  # param 1: alpha
  # param 2: beta
  # param 3: ratio gamma/sigma

  alpha <- param[1]
  beta  <- param[2]
  g_sigma <- param[3]

  logger::log_trace("Should go through tuning? {tuning}")
  logger::log_trace("Running for tune parameters: ",
                    "alpha: {alpha}, beta: {beta}, g_sigma: {g_sigma} ...")

  w_obs <- treatment_data

  # Normalization function with optional precomputed min and max
  normalize_min_max <- function(mat, precomputed_min = NULL, precomputed_max = NULL) {
    current_min <- ifelse(is.null(precomputed_min), min(mat), precomputed_min)
    current_max <- ifelse(is.null(precomputed_max), max(mat), precomputed_max)

    (mat - current_min) / (current_max - current_min)
  }

  # Now create distance matrix from spatial coordinates, scaled by beta
  loc_matrix <- as.matrix(stats::dist(spatial_coords))
  norm_loc_matrix <- normalize_min_max(loc_matrix)
  scaled_loc_matrix <- norm_loc_matrix * sqrt(1 / beta)

  # Store the min and max values for later use
  loc_min <- min(loc_matrix)
  loc_max <- max(loc_matrix)

  # Similarly for gps_matrix
  gps_matrix <- as.matrix(stats::dist(gps))
  norm_gps_matrix <- normalize_min_max(gps_matrix)
  scaled_gps_matrix <- norm_gps_matrix * sqrt(1 / alpha)

  gps_min <- min(gps_matrix)
  gps_max <- max(gps_matrix)

  # Create min_max table
  min_max_table <- data.frame(var = c("loc", "gps"), min = c(loc_min, gps_min),
                              max = c(loc_max, gps_max))


  # Here the scaled observations are not actually scaled
  scaled_obs <- cbind(w_obs, gps)
  colnames(scaled_obs) <- c("w_obs", "gps_obs")


  # Combine the squared matrices and then take the square root to get the Euclidean distance
  # between spatial location and gps
  combined_dist_matrix <- sqrt(scaled_loc_matrix^2 + scaled_gps_matrix^2)

  # Why is the "nugget" term set to be 1 here? Should it be smaller, or should this be tuned?
  t_sigma_obs_1 <- proc.time()
  sigma_obs <- g_sigma * kernel_fn(combined_dist_matrix) +
    diag(nrow(combined_dist_matrix))
  t_sigma_obs_2 <- proc.time()

  logger::log_trace("Wall clock time to generate covariance matrix ",
                    "({nrow(sigma_obs)},{ncol(sigma_obs)}): ",
                    "{t_sigma_obs_2[[3]] - t_sigma_obs_1[[3]]} s.")

  inv_sigma_obs <- compute_inverse(sigma_obs)

  # Estimate noise
  # Curious how this is used in the GP framework
  if (!tuning) {
    noise_est <- estimate_noise_gp(data = outcome_data,
                                   sigma_obs = sigma_obs,
                                   inv_sigma_obs = inv_sigma_obs)
    logger::log_debug("Estimated noise: {noise_est} ")
  }

  logger::log_debug("Computing weight and covariate balance for each ",
                    "requested exposure value ... ")

  col_all_list <- lapply(w,
                         function(w_instance) {

                           # compute weights
                           weights_res <- compute_weight_gp_spatial(w = w_instance,
                                                                    w_obs = w_obs,
                                                                    scaled_obs = scaled_obs,
                                                                    spatial_coords = spatial_coords,
                                                                    min_max_table = min_max_table,
                                                                    hyperparam = hyperparam,
                                                                    inv_sigma_obs = inv_sigma_obs,
                                                                    gps_m = gps_m,
                                                                    est_sd = !tuning,
                                                                    kernel_fn = kernel_fn)

                           weights_final <- weights_res$weight
                           weights_final[weights_final < 0] <- 0
                           if (sum(weights_final) > 0) {
                             weights_final <- weights_final / sum(weights_final)
                           }

                           # weigts.final = invers of paranthesis * kappa
                           # est is the same as m in the paper.

                           if (!tuning) {
                             est <- outcome_data %*% weights_final
                             pst_sd <- noise_est * weights_res$sd_scaled
                             logger::log_trace("Posterior for w = {w_instance} ==> ",
                                               "mu: {est}, var:{pst_sd}")
                           } else {
                             est <- NA
                             pst_sd <- NA
                           }
                           cov_balance_obj <- compute_w_corr(w = treatment_data,
                                                             covariate = covariates_data,
                                                             weight = weights_final)
                           covariate_balance <- as.vector(cov_balance_obj$absolute_corr)
                           c(covariate_balance, est, pst_sd)
                           list(covariate_balance = cov_balance_obj,
                                est = est,
                                pst_sd = pst_sd)
                         })

  logger::log_debug("Done with computing weight and covariate balance for ",
                    "each requested exposure value. ")

  col_all <- sapply(col_all_list, function(x) {
    x$covariate_balance$absolute_corr
  })
  est <- sapply(col_all_list, function(x) x$est)
  pst_sd <- sapply(col_all_list, function(x) x$pst_sd)

  # compute original covariate balance of data
  if (!tuning) {
    cov_balance_obj_org <- compute_w_corr(w = treatment_data,
                                          covariate = covariates_data,
                                          weight = rep(1, nrow(covariates_data)))
    cb_org <- cov_balance_obj_org$absolute_corr
  } else {
    cb_org <- NA
  }

  if (!is.matrix(col_all)){
    # in case of one covariate col_all returns vector instead of matrix
    row_name <- names(col_all)[1]
    col_all <- matrix(col_all, nrow = 1)
    rownames(col_all) <- row_name
  }


  if (nrow(col_all) == 1){
    col_all_w_average <- mean(col_all, na.rm = TRUE)
    names(col_all_w_average) <- rownames(col_all)
  } else {
    col_all_w_average <- rowMeans(col_all, na.rm = TRUE)
  }

  list(cb = col_all_w_average,
       cb_org = cb_org,
       est = est,
       pst = pst_sd)
}
