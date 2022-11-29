test_that("estimate_noise_gp works as expected", {

  set.seed(1099)
  data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)


  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  GPS <- GPS_m$GPS
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std

  kernel_fn <- function(x) exp(-x^2)
  hyperparam <- c(0.1, 0.2, 1)

  alpha <- hyperparam[1]
  beta  <- hyperparam[2]
  g_sigma <- hyperparam[3]

  w_obs <- data[[2]]

  scaled_obs <- cbind(w_obs*sqrt(1/alpha), GPS*sqrt(1/beta))
  sigma_obs <- g_sigma*kernel_fn(as.matrix(dist(scaled_obs))) + diag(nrow(scaled_obs))
  inv_sigma_obs <- compute_inverse(sigma_obs)

  noise_est <- estimate_noise_gp(data = data,
                                 sigma_obs = sigma_obs,
                                 inv_sigma_obs = inv_sigma_obs)

  expect_equal(length(noise_est), 1)
  #expect_equal(noise_est, 29.34867, tolerance = 0.0001)
})
