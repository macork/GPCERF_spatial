test_that("find_optimal_nn works as expected!", {
  set.seed(89)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_GPS(cov_mt = as.matrix(data[,-(1:2)]),
                     w_all = as.matrix(data$treat))

  # Hyperparameter
  hyperparam <- c(0.1, 0.2, 1)
  n_neighbor <- 10
  expand <- 1
  block_size <- 10000

  # compute posterior mean and standard deviation for vector of w.
  w <- seq(0,20,2)
  design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])

  hyperparam_grid <- expand.grid(seq(0.5,2.5,1),
                                 seq(0.4,0.6,0.2),
                                 seq(0.5,1.5,1))

  optimal_cb <- find_optimal_nn(w_obs = data$treat,
                                w = w,
                                y_obs = data$Y,
                                GPS_m = GPS_m,
                                design_mt = design_mt,
                                hyperparams = hyperparam_grid,
                                n_neighbor = 50, expand = 2, block_size = 2e3)

  opt_idx_nn <- order(colMeans(abs(optimal_cb)))[1]
  nn_opt_param <- unlist(hyperparam_grid[opt_idx_nn,])

  expect_equal(nrow(optimal_cb), 6L)
  expect_equal(ncol(optimal_cb), 12L)
  expect_equal(sum(is.na(optimal_cb)), 0)

  expect_equal(nn_opt_param[[1]], 0.5)
  expect_equal(nn_opt_param[[2]], 0.6)
  expect_equal(nn_opt_param[[3]], 0.5)
})
