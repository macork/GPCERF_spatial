test_that("compute_deriv_nn works as expected!", {

  set.seed(765)
  data <- generate_synthetic_data(sample_size = 200)
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  wi <- 4.8

  deriv_val <- compute_deriv_nn(w = wi,
                                w_obs = data$treat,
                                GPS_m = GPS_m,
                                y_obs = data$Y,
                                hyperparam = c(0.1,0.2,1),
                                n_neighbor = 20,
                                expand = 1,
                                block_size = 1000)

  expect_equal(length(deriv_val), 1L)
  expect_equal(deriv_val[1,1], -0.1793515, tolerance = 0.00001)
  expect_true(is.matrix(deriv_val))

})
