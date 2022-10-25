#' @title
#' Calculate Right Minus Left Derivatives for Change-point Detection in nnGP
#'
#' @description
#' Calculates the posterior mean of the difference between left- and
#' right-derivatives at an exposure level for the detection of change points.
#' nnGP approximation is used.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS values.
#'   - Column 2: Prediction of exposure for covariate of each data
#'   sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std).
#' @param y_obs A vector of observed outcome values.
#' @param hyperparam A vector of hyper-parameters in the GP model.
#' @param n_neighbor The number of nearest neighbors on one side
#' (see also \code{expand}).
#' @param expand A scaling factor to determine the total number of nearest
#' neighbors. The total is \code{2*expand*n_neighbor}.
#' @param block_size The number of samples included in a computation block.
#' Mainly used to balance the speed and memory requirement. Larger
#' \code{block_size} is faster, but requires more memory.
#' @param kernel_fn The covariance function. The input is the square of
#' Euclidean distance.
#' @param kernel_deriv_fn The partial derivative of the covariance function.
#' The input is the square of Euclidean distance.
#'
#' @return
#' A numeric value of the posterior mean of the difference between two one-sided
#' derivatives.
#'
#' @export
#'
#' @examples
#'
#' set.seed(325)
#' data <- generate_synthetic_data(sample_size = 200)
#' GPS_m <- train_GPS(cov_mt = data[,-(1:2)], w_all = data$treat)
#'
#' wi <- 12.2
#'
#' deriv_val <- compute_rl_deriv_nn(w = wi,
#'                                  w_obs = data$treat,
#'                                  GPS_m = GPS_m,
#'                                  y_obs = data$Y,
#'                                  hyperparam = c(0.2,0.4,1.2),
#'                                  n_neighbor = 20,
#'                                  expand = 1,
#'                                  block_size = 1000)
#'
compute_rl_deriv_nn <-  function(w,
                                 w_obs,
                                 GPS_m,
                                 y_obs,
                                 hyperparam,
                                 n_neighbor,
                                 expand,
                                 block_size,
                                 kernel_fn = function(x) exp(-x),
                                 kernel_deriv_fn = function(x) -exp(-x)
                                 ){

  left_deriv <- compute_deriv_nn(w,
                                 w_obs[w_obs<w],
                                 GPS_m[w_obs<w,],
                                 y_obs[w_obs<w],
                                 hyperparam,
                                 n_neighbor = n_neighbor,
                                 expand = expand,
                                 block_size = block_size,
                                 kernel_fn = kernel_fn,
                                 kernel_deriv_fn = kernel_deriv_fn)

  right_deriv <- compute_deriv_nn(w,
                                  w_obs[w_obs>=w],
                                  GPS_m[w_obs>=w,],
                                  y_obs[w_obs>=w],
                                  hyperparam,
                                  n_neighbor = n_neighbor,
                                  expand = expand,
                                  block_size = block_size,
                                  kernel_fn = kernel_fn,
                                  kernel_deriv_fn = kernel_deriv_fn)

  return(right_deriv - left_deriv)
}
