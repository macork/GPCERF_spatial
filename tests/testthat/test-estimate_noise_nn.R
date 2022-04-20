test_that("estimate_noise_nn works as expected!", {

  set.seed(425)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                     w.all = as.matrix(data$treat))

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

  noise <- estimate_noise_nn(hyperparam = hyperparam,
                             w_obs = data$treat,
                             GPS_obs = GPS_m$GPS,
                             y_obs = y_use_ord,
                             n_neighbor = n_neighbor)


  expect_equal(noise, 31.34945, tolerance = 0.0001)
})
