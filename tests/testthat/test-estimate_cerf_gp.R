test_that("estimate_cerf_gp works as expected!", {

  set.seed(129)
  sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  # In the future, CausalGPS gps estimation will be used.
  GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                     w.all = as.matrix(sim.data$treat))

  # exposure values
  w.all <- seq(0,20,0.1)

  data.table::setDT(sim.data)

  cerf_gp_obj <- estimate_cerf_gp(data = sim.data,
                                  w = w.all,
                                  GPS_m = GPS_m,
                                  params = list(alpha = c(0.1,0.2,0.4),
                                                beta=0.2,
                                                g_sigma = 1,
                                                tune_app = "all"),
                                  nthread = 1)


  # estimate_cerf_gp returns S3 class cerf_gp
  expect_s3_class(cerf_gp_obj, "cerf_gp")

  expect_equal(length(cerf_gp_obj$pst_mean), 201L)
  expect_equal(length(cerf_gp_obj$w), 201L)
  expect_equal(cerf_gp_obj$pst_mean[1], -13.54376, tolerance = 0.00001)
  expect_equal(cerf_gp_obj$pst_mean[10], -25.41547, tolerance = 0.00001)
  expect_equal(cerf_gp_obj$w[70], w.all[70], tolerance = 0.00001)
})
