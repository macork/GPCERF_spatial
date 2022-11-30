set.seed(86)
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
# Estimate GPS function
GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                  w_all = data$treat,
                  sl_lib = c("SL.xgboost"),
                  dnorm_log = FALSE)
# Hyperparameter
hyperparam <- c(0.1, 0.2, 1)
n_neighbor <- 15
expand <- 1
block_size <- 10000
# compute noise
noise <- estimate_noise_nn(hyperparam = hyperparam,
                          w_obs = data$treat,
                          GPS_obs = GPS_m$GPS,
                          y_obs = data$Y,
                          n_neighbor = n_neighbor)
# compute posterior mean and standard deviation for vector of w.
w <- seq(0,20,1)
val <- estimate_mean_sd_nn(hyperparam = hyperparam,
                          sigma2 = noise,
                          w_obs = data$treat,
                          w = w,
                          y_obs = data$Y,
                          GPS_m = GPS_m,
                          n_neighbor = n_neighbor,
                          expand = expand,
                          block_size = block_size,
                          nthread = 1)
