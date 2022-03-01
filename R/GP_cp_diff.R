#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param w
#' @param y.obs
#' @param w.obs
#' @param GPS.obs
#' @param param
#' @param e_gps_pred
#' @param e_gps_std
#' @param noise_est
#'
#' @return
#' @export
#'
#' @examples
GP.cp.diff = function(w, y.obs, w.obs, GPS.obs, param, e_gps_pred, e_gps_std, noise_est){
  # left derivatives
  left.weights = GP.deriv.weights.calc(w, w.obs[w.obs<w], GPS.obs[w.obs<w],
                                       param, e_gps_pred, e_gps_std)
  right.weights = GP.deriv.weights.calc(w, w.obs[w.obs>=w], GPS.obs[w.obs>=w],
                                        param, e_gps_pred, e_gps_std)
  right.weights%*%y.obs[w.obs>=w] - left.weights%*%y.obs[w.obs>=w]
}
