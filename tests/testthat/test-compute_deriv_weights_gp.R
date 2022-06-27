test_that("compute_deriv_weights_gp works as expected!", {

  set.seed(9615)
  data <- generate_synthetic_data(sample_size = 200)
  GPS_m <- train_GPS(cov_mt = as.matrix(data[,-(1:2)]),
                     w_all = as.matrix(data$treat))

  wi <- 4.2
  weights <- compute_deriv_weights_gp(w = wi,
                                      w_obs = data$treat,
                                      GPS_m = GPS_m,
                                      hyperparam = c(1,1,2))

  expect_equal(length(weights), nrow(data))
  expect_equal(weights[37], -3.169395e-04, tolerance = 0.00001)
})
