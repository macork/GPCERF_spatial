test_that("compute_posterior_sd_nn works as expected.", {

  set.seed(959)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_GPS(cov_mt = data[,-(1:2)],
                     w_all = data$treat)

  # Hyperparameter
  hyperparam <- c(0.1, 0.2, 1)
  n_neighbor <- 10
  expand <- 1
  block_size <- 10000

  # Exposure level
  wi <- 0.4

  # Estimate GPS for the exposure level
  GPS_w = dnorm(wi,
                mean = GPS_m$e_gps_pred,
                sd = GPS_m$e_gps_std, log = TRUE)

  # Order data for easy selection
  coord_obs = cbind(data$treat, GPS_m$GPS)
  y_use <- data$Y

  obs_ord <- coord_obs[order(coord_obs[,1]),]
  y_use_ord <- y_use[order(coord_obs[,1])]

  # compute noise
  noise <- estimate_noise_nn(hyperparam = hyperparam,
                             w_obs = data$treat,
                             GPS_obs = GPS_m$GPS,
                             y_obs = y_use_ord,
                             n_neighbor = n_neighbor)

  # compute posterior standard deviation
  pst_sd <- compute_posterior_sd_nn(hyperparam = hyperparam,
                                    w = wi,
                                    GPS_w = GPS_w,
                                    obs_ord = obs_ord,
                                    sigma2 = noise,
                                    n_neighbor = 20,
                                    expand = 1)

  # expect_equal(pst_sd, 5.437376, tolerance = 0.00001)
})
