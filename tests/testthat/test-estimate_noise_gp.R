test_that("estimate_noise_gp works as expected", {

  set.seed(1099)
  data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)
  data.table::setDT(data)

  # Estimate GPS function
  GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                     w.all = as.matrix(data$treat))

  hyperparam <- c(0.1, 0.2, 1)
  noise_est <- estimate_noise_gp(hyperparam, data, GPS_m$GPS)

  expect_equal(length(noise_est), 1)
  expect_equal(noise_est, 29.34867, tolerance = 0.0001)
})
