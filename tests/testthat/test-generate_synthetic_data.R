test_that("generate_synthetic_data works as expected!", {

  set.seed(298)
  mydata <- generate_synthetic_data(sample_size = 200)

  expect_equal(class(mydata), "data.frame")
  expect_equal(nrow(mydata), 200L)
  expect_equal(length(mydata), 8L)
  expect_equal(mydata[10,2], 9.771912, tolerance = 0.00001)

})
