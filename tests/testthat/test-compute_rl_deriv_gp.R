test_that("compute_rl_deriv_gp works as expected!", {

  set.seed(127)
  data <- generate_synthetic_data(sample_size = 200)
  GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                     w.all = as.matrix(data$treat))

  wi <- 8.6

  deriv_val <- compute_rl_deriv_gp(w = wi,
                                   w_obs = data$treat,
                                   y_obs = data$Y,
                                   GPS_m = GPS_m,
                                   hyperparam = c(0.2,0.8,2))

  expect_equal(length(deriv_val), 1L)
  expect_equal(deriv_val[1,1], -0.003027017, tolerance = 0.00001)
  expect_true(is.matrix(deriv_val))

})
