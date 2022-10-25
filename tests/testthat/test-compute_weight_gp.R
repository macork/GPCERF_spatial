test_that("multiplication works", {

  set.seed(917)
  # Generate synthetic data
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
  w_obs <- obs_exposure <- data$treat

  # Choose an exposure level to compute CERF
  w = 1.8

  # Define kernel function
  kernel_fn = function(x) exp(-x^2)

  # compute GPS, e_gps_pred, and e_gps_std
  e_gps <- xgboost(label=data$treat, data=as.matrix(data[,-(1:2)]),
                   nrounds = 50)
  e_gps_pred <- predict(e_gps,as.matrix(data[,-(1:2)]))
  e_gps_std <- sd(data$treat-e_gps_pred)
  GPS <- dnorm(data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)
  GPS_m <- data.frame(GPS, e_gps_pred, e_gps_std)

  # set hyperparameters
  hyperparam <- c(0.1, 0.4, 1)
  alpha <- hyperparam[1]
  beta <- hyperparam[2]
  g_sigma <- hyperparam[3]

  # Compute scaled observation data and inverse of covariate matrix.
  scaled_obs = cbind(obs_exposure*sqrt(1/alpha), GPS*sqrt(1/beta))
  sigma_obs = g_sigma*kernel_fn(as.matrix(dist(scaled_obs))) + diag(nrow(scaled_obs))
  inv_sigma_obs <- compute_inverse(sigma_obs)


  weight <- compute_weight_gp(w = w,
                              w_obs = w_obs,
                              scaled_obs = scaled_obs,
                              hyperparam = hyperparam,
                              inv_sigma_obs = inv_sigma_obs,
                              GPS_m = GPS_m,
                              kernel_fn = kernel_fn)

  expect_equal(length(weight), 200L)
  expect_equal(weight[28], 1.588098e-03, tolerance = 10e-5)
})
