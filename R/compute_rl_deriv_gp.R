#' @title
#' Change-point Detection in Full GP
#'
#' @description
#' Calculates the posterior mean of the difference between left- and
#' right-derivatives at an exposure level for the detection of change points.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param y_obs A vector of observed outcome values of all samples.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param hyperparam A vector of hyper-parameters in the GP model.
#' @param kernel_fn The covariance function.
#' @param kernel_deriv_fn The partial derivative of the covariance function.
#'
#' @return
#' A numeric value of the posterior mean of the difference between two one-sided
#' derivatives.
#' @export
#'
#' @examples
#'
#' set.seed(847)
#' data <- generate_synthetic_data(sample_size = 200)
#' GPS_m <- train_GPS(cov_mt = data[,-(1:2)], w_all = data$treat)
#'
#' wi <- 8.6
#'
#' val <- compute_rl_deriv_gp(w = wi,
#'                            w_obs = data$treat,
#'                            y_obs = data$Y,
#'                            GPS_m = GPS_m,
#'                            hyperparam = c(1,1,2))
#'
compute_rl_deriv_gp <- function(w,
                                w_obs,
                                y_obs,
                                GPS_m,
                                hyperparam,
                                kernel_fn = function(x) exp(-x),
                                kernel_deriv_fn = function(x) -exp(-x)){
  # left side weights
  left_weights <-  compute_deriv_weights_gp(w = w,
                                            w_obs = w_obs[w_obs<w],
                                            GPS_m = GPS_m[w_obs<w,],
                                            hyperparam = hyperparam,
                                            kernel_fn = kernel_fn,
                                            kernel_deriv_fn = kernel_deriv_fn)

  # right side weights
  right_weights <-  compute_deriv_weights_gp(w = w,
                                             w_obs = w_obs[w_obs>=w],
                                             GPS_m = GPS_m[w_obs>=w,],
                                             hyperparam = hyperparam,
                                             kernel_fn = kernel_fn,
                                             kernel_deriv_fn = kernel_deriv_fn)

  # compute derivative
  return(right_weights%*%y_obs[w_obs>=w] - left_weights%*%y_obs[w_obs<w])
}
