#' @title
#' Calculates Right Minus Left Derivatives for Change-point Detection in nnGP
#'
#' @description
#' Calculates the posterior mean of the difference between left- and right-derivatives
#' at an exposure level for the detection of change points. nnGP approximation is used.
#'
#' @param w A numeric value of the exposure level at which the difference between two one-sided
#' partial derivatives is calculated.
#' @param w.obs A vector of observed exposure levels.
#' @param GPS.obs A vector of estimated GPS at the observed exposure levels.
#' @param y.obs A vector of observed outcomes.
#' @param param A vector of hyperparameters for the GP.
#' @param e_gps_pred A vector of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param kernel.fn The covariance function of GP.
#' @param kernel.deriv.fn The partial derivative of the covariance function.
#' @param n.neighbour Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n.neighbour}.
#'
#' @return A numeric value of the posterior mean of the difference between two one-sided
#' derivatives.
#' @export
#'
#' @examples
nn.cp.calc = function(w, w.obs, GPS.obs, y.obs, param, e_gps_pred, e_gps_std,
                      kernel.fn = function(x) exp(-x^2),
                      kernel.deriv.fn = function(x, mu, sigma) (x-mu)/sigma^2,
                      n.neighbour, expand, block.size){
  left.deriv = deriv.nn.fast(w, w.obs[w.obs<w], GPS.obs[w.obs<w], y.obs[w.obs<w], param,
                             e_gps_pred, e_gps_std, n.neighbour = n.neighbour,
                             expand = expand, block.size = block.size,
                             kernel.fn = kernel.fn, kernel.deriv.fn = kernel.deriv.fn)
  right.deriv = deriv.nn.fast(w, w.obs[w.obs>=w], GPS.obs[w.obs>=w], y.obs[w.obs>=w], param,
                               e_gps_pred, e_gps_std, n.neighbour = n.neighbour,
                               expand = expand, block.size = block.size,
                               kernel.fn = kernel.fn, kernel.deriv.fn = kernel.deriv.fn)

  right.deriv - left.deriv
}
