#' @title
#' Calculate weights for estimation of a point on CERF
#'
#' @description
#' Calculates the weights of observed outcomes which is then used to estimate
#' the posterior mean of CERF at a given exposure level.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param scaled_obs A matrix of two columns.
#'   - First column is the scaled GPS value of all samples
#'   (GPS * 1 / sqrt(alpha))
#'   - Second column is the scaled exposure value of all samples
#'   (w * 1/sqrt(beta))
#' @param hyperparam A vector of hyper-parameters for the GP.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: gamma/sigma
#' @param inv_sigma_obs Inverse of the covariance matrix between observed
#' samples.
#' @param gps_m An S3 gps object including:
#'   gps: A data.frame of GPS vectors.
#'     - Column 1: GPS
#'     - Column 2: Prediction of exposure for covariate of each data sample
#'     (e_gps_pred).
#'     - Column 3: Standard deviation of  e_gps (e_gps_std)
#'   used_params:
#'     - dnorm_log: TRUE or FLASE
#' @param est_sd Should the posterior se be computed (default=FALSE)
#' @param kernel_fn The covariance function of GP.
#'
#' @return
#' A list of two elements, weight and standard deviation.
#'
#' @keywords internal
#'
compute_weight_gp_spatial <- function(w, w_obs, scaled_obs, spatial_coords, min_max_table,
                              hyperparam, inv_sigma_obs, gps_m, est_sd = FALSE,
                              kernel_fn = function(x) exp(-x ^ 2)) {

  logger::log_trace("Computing weights for w = {w} ...")

  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]

  # Compute GPS for requested w
  e_gps_pred <- gps_m$gps$e_gps_pred
  e_gps_std <- gps_m$gps$e_gps_std
  dnorm_log <- gps_m$used_params$dnorm_log
  gps_w <- stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = dnorm_log)

  # Extract normalization parameters
  loc_min <- min_max_table$min[min_max_table$var == "loc"]
  loc_max <- min_max_table$max[min_max_table$var == "loc"]
  gps_min <- min_max_table$min[min_max_table$var == "gps"]
  gps_max <- min_max_table$max[min_max_table$var == "gps"]

  # Filtering observed data based on the exposure constraint
  # Current constraint is to be within 3 of actual observed
  valid_obs_indices <- which(abs(w_obs - w) <= 4)
  filtered_scaled_obs <- scaled_obs[valid_obs_indices, ]
  filtered_spatial_coords <- spatial_coords[valid_obs_indices, ]

  n <- length(gps_w)

  gps_matrix <- matrix(0, n, n)
  loc_matrix <- matrix(0, n, n)

  # Create gps and scaled loc matrix excluding those farther than 4 away.
  # Replace this with more efficient process eventually
  for(i in 1:n) {
    for(j in 1:n) {
      # For GPS
      if(j %in% valid_obs_indices) {
        gps_matrix[i,j] <- abs(gps_w[i] - scaled_obs[j, "gps_obs"])
      } else {
        gps_matrix[i,j] <- Inf  # Set to a large value
      }

      # For spatial coords
      if(j %in% valid_obs_indices) {
        loc_matrix[i,j] <- sqrt(sum((spatial_coords[i,] - spatial_coords[j,])^2))
      } else {
        loc_matrix[i,j] <- Inf  # Set to a large value
      }
    }
  }

  # Normalize using computed min/max
  norm_gps_matrix <- normalize_min_max(gps_matrix, gps_min, gps_max)

  # Scale by appropriate ammount
  scaled_gps_matrix <- norm_gps_matrix * sqrt(1 / alpha)

  # Normalize using computed min/max
  norm_loc_matrix <- normalize_min_max(loc_matrix, loc_min, loc_max)

  # Scale by appropriate amount
  scaled_loc_matrix <- norm_loc_matrix * sqrt(1 / beta)

  # Now combine
  combined_dist_matrix <- sqrt(scaled_loc_matrix^2 + scaled_gps_matrix^2)

  # kappa
  # sigma_cross = kappa/sigma^2 : Is always n*n matrix.
  # each column of sigma_cross is ki.
  # statspat.geom::crossdist

  # kappa and sigma_cross

  sigma_cross <- g_sigma * kernel_fn(combined_dist_matrix)

  logger::log_trace("sigma_cross {class(sigma_cross)[1]} ",
                    "({nrow(sigma_cross)}, {ncol(sigma_cross)}) ",
                    "was generated.")

  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  normalized_sigma_cross <- Rfast::colmeans(sigma_cross)
  weight <- c(arma_mm(inv_sigma_obs, normalized_sigma_cross))

  # Make weights 0 for those too far away from cutoff (usually negative weights)
  weight[-valid_obs_indices] <- 0

  logger::log_trace("Weight vector with size {length(weight)}",
                    " was generated.")

  # compute scaled posterior sd
  if (est_sd) {

    # Create GPS diff with itself
    gps_diff <- outer(gps_w, gps_w, "-")
    # Replace with Inf for invalid indices
    gps_diff[, -valid_obs_indices] <- Inf
    gps_self_matrix <- abs(gps_diff)

    # Normalize using computed min/max
    norm_gps_self_matrix <- normalize_min_max(gps_self_matrix, gps_min, gps_max)

    # Scale by appropriate ammount
    scaled_gps_self_matrix <- norm_gps_matrix * sqrt(1 / alpha)

    # Now combine
    combined_dist_matrix <- sqrt(scaled_loc_matrix^2 + scaled_gps_self_matrix^2)

    # This not correct, but using for now
    sigma_w <- g_sigma * kernel_fn(combined_dist_matrix) +
      diag(n)
    sd_scaled <- sqrt(sum(sigma_w) / n ^ 2 -
                        sum(weight * normalized_sigma_cross))
    sd_scaled = 1
    logger::log_trace("Computed scaled standard deviation: {sd_scaled}")
  } else {
    sd_scaled <- NA
  }

  return(list(weight = weight, sd_scaled = sd_scaled))
}
