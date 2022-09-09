#' @title
#' Calculate Weights for Estimation of a Point on CERF
#'
#' @description
#' Calculates the weights of observed outcomes which is then used to estimate
#' the posterior mean of CERF at a given exposure level.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param scaled_obs A matrix of two columns.
#'   - First column is the scaled GPS value of all samples (GPS * 1/sqrt(alpha))
#'   - Second column is the scaled exposure value of all samples (w * 1/sqrt(beta))
#' @param hyperparam A vector of hyper-parameters for the GP.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: gamma/sigma
#' @param inv_sigma_obs Inverse of the covariance matrix between observed samples.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: A vector of estimated GPS evaluated at the observed exposure levels.
#'   - Column 2: Estimated conditional means of the exposure given covariates
#'               for all samples (e_gps_pred).
#'   - Column 3: Estimated conditional standard deviation of the exposure given
#'               covariates for all samples (e_gps_std).
#' @param kernel_fn The covariance function of GP.
#'
#' @return
#' A vector of the weights assigned to each sample for the calculate of posterior
#'  mean of CERF at \code{w}.
#' @export
#'
#' @examples
#'
#' set.seed(814)
#' #Generate synthetic data
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#' w_obs <- obs_exposure <- data$treat
#'
#' # Choose an exposure level to compute CERF
#' w = 1.8
#'
#' # Define kernel function
#' kernel_fn <- function(x) exp(-x^2)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' GPS <- GPS_m$GPS
#'
#' # set hyperparameters
#' hyperparam <- c(0.1, 0.4, 1)
#' alpha <- hyperparam[1]
#' beta <- hyperparam[2]
#' g_sigma <- hyperparam[3]
#'
#' # Compute scaled observation data and inverse of covariate matrix.
#' scaled_obs <- cbind(obs_exposure*sqrt(1/alpha), GPS*sqrt(1/beta))
#' sigma_obs <- g_sigma*kernel_fn(as.matrix(dist(scaled_obs))) + diag(nrow(scaled_obs))
#' inv_sigma_obs <- compute_inverse(sigma_obs)
#'
#'
#' weight <- compute_weight_gp(w = w,
#'                             w_obs = w_obs,
#'                             scaled_obs = scaled_obs,
#'                             hyperparam = hyperparam,
#'                             inv_sigma_obs = inv_sigma_obs,
#'                             GPS_m = GPS_m,
#'                             kernel_fn = kernel_fn)
#'
#'
compute_weight_gp <- function(w, w_obs, scaled_obs, hyperparam,
                              inv_sigma_obs, GPS_m,
                              kernel_fn = function(x) exp(-x^2)){

  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]

  # Compute GPS for requested w
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std

  # TODO: The following section is repeated between this function
  # and compute_sd_gp function.
  GPS_w <- stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)
  scaled_w <- cbind( w*sqrt(1/alpha), GPS_w*sqrt(1/beta) )

  # kappa
  # sigma_cross = kappa/sigma^2 : Is always n*n matrix.
  # each column of sigma_cross is ki.
  # statspat.geom::crossdist
  sigma_cross <- g_sigma*kernel_fn(crossdist(scaled_w[,1],
                                             scaled_w[,2],
                                             scaled_obs[,1],
                                             scaled_obs[,2]))

  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  normalized_sigma_cross <- Rfast::colmeans(sigma_cross)   #rep(1/length(w_obs),length(w_obs))%*%sigma_cross
  weight <- c(arma_mm(inv_sigma_obs, normalized_sigma_cross))   #c((normalized_sigma_cross)%*%inv_sigma_obs)

  return(weight)
}
