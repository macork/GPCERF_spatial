test_that("calc_ac works as expected!", {

  set.seed(729)

  # generate data
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # generate random weights
  weights <- runif(nrow(data))
  weights <- weights / sum(weights)

  # covariate matrix
  design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])

  cb <- calc_ac(w = data$treat, X = design_mt, weights=weights)

  expect_equal(length(cb), 6L)
  expect_equal(cb[1], 0.06868856, tolerance = 0.00001)
  expect_equal(cb[5], 0.31722877, tolerance = 0.00001)
})
