test_that("estimate_cerf_nngp works as expected!", {

  set.seed(19)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
  # Estimate GPS function
  GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                     w_all = data$treat,
                     sl_lib = c("SL.xgboost"),
                     dnorm_log = FALSE)

  # exposure values
  w_all <- seq(0,20,0.5)
  cerf_nngp_obj <- estimate_cerf_nngp(data,
                                      w_all,
                                      GPS_m,
                                      params = list(alpha = c(0.1,0.2),
                                                    beta = 0.2,
                                                    g_sigma = 1,
                                                    tune_app = "all",
                                                    n_neighbor = 20,
                                                    expand = 1,
                                                    block_size = 1e4),
                                      formula = ~ . - 1 - Y - treat)

  expect_s3_class(cerf_nngp_obj, "cerf_nngp")

  expect_equal(length(cerf_nngp_obj$pst_mean), 41L)
  expect_equal(length(cerf_nngp_obj$w), 41L)
  #expect_equal(cerf_nngp_obj$pst_mean[1], -22.21065, tolerance = 0.00001)
  #expect_equal(cerf_nngp_obj$pst_mean[10], -8.955434, tolerance = 0.00001)
  expect_equal(cerf_nngp_obj$w[31], w_all[31], tolerance = 0.00001)

})
