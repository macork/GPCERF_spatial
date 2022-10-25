test_that("Weighted correlation works as expected.", {

  set.seed(239)
  data1 <- generate_synthetic_data(sample_size = 100)
  weights1 <- runif(nrow(data1))
  val1 <- compute_w_corr(data = data1, weights = weights1)

  expect_vector(val1)
  expect_equal(length(val1), 6L)
  expect_equal(val1[1], 0.9725398, tolerance = 0.00001)
  expect_equal(val1[2], 0.0259583, tolerance = 0.00001)


  # number of data.samples and weights should be the same
  data3 <- generate_synthetic_data(sample_size = 50)
  weights3 <- runif(nrow(data3)+20)
  expect_error(compute_w_corr(data = data3, weights = weights3))

})
