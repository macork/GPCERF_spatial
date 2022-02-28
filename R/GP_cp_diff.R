GP.cp.diff = function(w, y.obs, w.obs, GPS.obs, param, e_gps_pred, e_gps_std, noise_est){
  # left derivatives
  left.weights = GP.deriv.weights.calc(w, w.obs[w.obs<w], GPS.obs[w.obs<w],
                                       param, e_gps_pred, e_gps_std)
  right.weights = GP.deriv.weights.calc(w, w.obs[w.obs>=w], GPS.obs[w.obs>=w],
                                        param, e_gps_pred, e_gps_std)
  right.weights%*%y.obs[w.obs>=w] - left.weights%*%y.obs[w.obs>=w]
}
