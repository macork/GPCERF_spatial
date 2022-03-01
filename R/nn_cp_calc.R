nn.cp.calc = function(w, GPS.new, y.obs, w.obs, GPS.obs, param, e_gps_pred, e_gps_std,
                      kernel.fn = function(x) exp(-x^2),
                      kernel.deriv.fn = function(x, mu, sigma) (x-mu)/sigma^2){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  left.deriv = deriv.nn.fast(params, w, GPS.new, cbind(w.obs,GPS.obs)[w.obs<w,], y.obs, n.neighbour = 200)
  right.deriv = deriv.nn.fast(params, w, GPS.new, cbind(w.obs,GPS.obs)[w.obs>=w,], y.obs, n.neighbour = 200)
  right.deriv - left.deriv
}
