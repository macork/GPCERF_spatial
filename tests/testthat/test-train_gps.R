test_that("train_gps works as expected.", {

  set.seed(651)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  expect_s3_class(GPS_m, "gps")
  expect_false(GPS_m$used_params$dnorm_log)
  expect_equal(GPS_m$gps[4, 1], 8.74310154, tolerance = 0.000001)
})
