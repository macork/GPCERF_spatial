test_that("estimate_cerf_gp works as expected!", {

  set.seed(129)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  # exposure values
  w_all <- seq(0, 20, 0.1)

  cerf_gp_obj <- estimate_cerf_gp(data = data,
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all"),
                                  nthread = 1)


  expect_s3_class(cerf_gp_obj, "cerf_gp")

  expect_equal(length(cerf_gp_obj$pst_mean), 201L)
  expect_equal(length(cerf_gp_obj$w), 201L)

  # Check input parameters -----------------------------------------------------
  set.seed(129)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  # exposure values
  w_all <- seq(0, 20, 0.1)

  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = as.matrix(data),
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all"),
                                  nthread = 1))

  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = data,
                                  w = w_all,
                                  GPS_m = as.matrix(GPS_m),
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all"),
                                  nthread = 1))


  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = data,
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                ggggg_sigma = 1,
                                                tune_app = "all"),
                                  nthread = 1))

  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = data,
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(),
                                                beta = c(),
                                                g_sigma = ,
                                                tune_app = "all"),
                                  nthread = 1))


  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = data,
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "at_random"),
                                  nthread = 1))

  expect_error(cerf_gp_obj <- estimate_cerf_gp(
                                  data = data,
                                  w = w_all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1, 0.2, 0.4),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "xyz"),
                                  nthread = 1))


  # Check with two threads -----------------------------------------------------
  set.seed(129)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  # exposure values
  w_all <- seq(0, 30, 0.1)

  cerf_gp_obj_2 <- estimate_cerf_gp(
                                   data = data,
                                   w = w_all,
                                   GPS_m = GPS_m,
                                   params = list(alpha = c(0.1, 0.2, 0.4),
                                                 beta = 0.2,
                                                 g_sigma = 1,
                                                 tune_app = "all"),
                                   nthread = 2)

  expect_s3_class(cerf_gp_obj_2, "cerf_gp")

  expect_equal(length(cerf_gp_obj_2$pst_mean), 301L)
  expect_equal(length(cerf_gp_obj_2$w), 301L)

})
