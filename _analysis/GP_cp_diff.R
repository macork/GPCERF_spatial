#' @title
#' Change-point Detection in Full GP
#'
#' @description
#' Calculate the posterior mean of the difference between left- and right-derivatives
#' at an exposure level for the detection of change points.
#'
#' @param w A numeric value of the exposure level at which the difference between two one-sided
#' partial derivatives is calculated.
#' @param y.obs A vector of observed outcomes.
#' @param w.obs A vector of observed exposure levels.
#' @param GPS.obs A vector of estimated GPS at the observed exposure levels.
#' @param param A vector of hyperparameters for the GP.
#' @param e_gps_pred A vector of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param kernel.fn The covariance function of GP.
#' @param kernel.deriv.fn The partial derivative of the covariance function.
#'
#' @return A numeric value of the posterior mean of the difference between two one-sided
#' derivatives.
#' @export
#'
#' @examples
GP.cp.diff = function(w, y.obs, w.obs, GPS.obs, param, e_gps_pred, e_gps_std, kernel.fn, kernel.deriv.fn){
  # left derivatives
  left.weights = GP.deriv.weights.calc(w, w.obs[w.obs<w], GPS.obs[w.obs<w],
                                       param, e_gps_pred, e_gps_std, kernel.fn,
                                       kernel.deriv.fn)
  right.weights = GP.deriv.weights.calc(w, w.obs[w.obs>=w], GPS.obs[w.obs>=w],
                                        param, e_gps_pred, e_gps_std, kernel.fn,
                                        kernel.deriv.fn)
  right.weights%*%y.obs[w.obs>=w] - left.weights%*%y.obs[w.obs>=w]
}
