test_that("estimate_cerf_nngp works as expected!", {

  set.seed(19)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- estimate_gps(cov_mt = data[, -(1:2)],
                        w_all = data$treat,
                        sl_lib = c("SL.xgboost"),
                        dnorm_log = FALSE)

  # exposure values
  w_all <- seq(0, 20, 0.5)
  cerf_nngp_obj <- estimate_cerf_nngp(data,
                                      w_all,
                                      GPS_m,
                                      params = list(alpha = c(0.1, 0.2),
                                                    beta = 0.2,
                                                    g_sigma = 1,
                                                    tune_app = "all",
                                                    n_neighbor = 20,
                                                    block_size = 1e4),
                                      formula = ~ . - 1 - Y - treat)

  expect_error(estimate_cerf_nngp(data,
                                  w_all,
                                  GPS_m,
                                  params = list(alpha = c(0.1, 0.2),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "at_random",
                                                n_neighbor = 20,
                                                block_size = 1e4),
                                  formula = ~ . - 1 - Y - treat))

  expect_error(estimate_cerf_nngp(data,
                                  w_all,
                                  GPS_m,
                                  params = list(alpha = c(0.1, 0.2),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "abc",
                                                n_neighbor = 20,
                                                block_size = 1e4),
                                  formula = ~ . - 1 - Y - treat))

  expect_s3_class(cerf_nngp_obj, "cerf_nngp")

  expect_equal(length(cerf_nngp_obj$pst_mean), 41L)
  expect_equal(length(cerf_nngp_obj$w), 41L)
  expect_equal(cerf_nngp_obj$w[31], w_all[31], tolerance = 0.00001)

  # expect error with missing data
  data_na <- data
  data_na$cf2[3] <- NA
  expect_error(estimate_cerf_nngp(data_na,
                                  w_all,
                                  GPS_m,
                                  params = list(alpha = c(0.1, 0.2),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all",
                                                n_neighbor = 20,
                                                block_size = 1e4),
                                  formula = ~ . - 1 - Y - treat))

  # Check non-consistent data and GPS object. ----------------------------------
  # Different size
  set.seed(659)
  data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)
  w_all <- seq(0, 20, 0.1)
  # Estimate GPS function
  GPS_m <- estimate_gps(cov_mt = data[,-(1:2)],
                        w_all = data$treat,
                        sl_lib = c("SL.xgboost"),
                        dnorm_log = FALSE)

  GPS_mm <- GPS_m
  GPS_mm$gps <- GPS_mm$gps[1:99, ]

  expect_error(estimate_cerf_nngp(data,
                                  w_all,
                                  GPS_mm,
                                  params = list(alpha = c(0.1, 0.2),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all",
                                                n_neighbor = 20,
                                                block_size = 1e4),
                                  formula = ~ . - 1 - Y - treat))

  # Same size but different exposure values
  GPS_mx <- GPS_m
  GPS_mx$gps$w[4] <- GPS_mm$gps$w[4] + 2.5

  expect_error(estimate_cerf_nngp(data,
                                  w_all,
                                  GPS_mx,
                                  params = list(alpha = c(0.1, 0.2),
                                                beta = 0.2,
                                                g_sigma = 1,
                                                tune_app = "all",
                                                n_neighbor = 20,
                                                block_size = 1e4),
                                  formula = ~ . - 1 - Y - treat))
})
