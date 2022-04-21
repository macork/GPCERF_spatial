
set.seed(9615)
data <- generate_synthetic_data(sample_size = 200)
GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                   w.all = as.matrix(data$treat))

wi <- 4.2

weights = compute_deriv_weights_gp(w = wi,
                                   w_obs = data$treat,
                                   GPS_m = GPS_m,
                                   hyperparam = c(1,1,2))

deriv.est.full = sapply(seq(2.5,17.5,0.1), function(w){
  weights = GP.deriv.weights.calc(w = w,
                                  w_obs =
                                  sim.data$treat, GPS.obs = GPS,
                                  param = c(1,1,2), e_gps_pred = e_gps_pred,
                                  e_gps_std = e_gps_std)
  weights%*%sim.data$Y
})
