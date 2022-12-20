set.seed(915)
data <- generate_synthetic_data(sample_size = 150)
GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                  w_all = data$treat,
                  sl_lib = c("SL.xgboost"),
                  dnorm_log = FALSE)
wi <- 4.8
weights <- compute_deriv_weights_gp(w = wi,
                                   w_obs = data$treat,
                                   GPS_m = GPS_m,
                                   hyperparam = c(1,1,2))
