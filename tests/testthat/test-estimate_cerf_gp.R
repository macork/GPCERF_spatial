test_that("estimate_cerf_gp works as expected!", {

  sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)

  # Compute true curve
  tru.curve <- sapply(seq(0,20,0.1), function(w){
    tru_R(w, sim.data)
  })

  # Estimate GPS function
  # In the future, CausalGPS gps estimation will be used.
  GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                     w.all = as.matrix(sim.data$treat))

  # exposure values
  w.all = seq(0,20,0.1)

  data.table::setDT(sim.data)

  cerf_gp_obj <- estimate_cerf_gp(sim.data,
                                  w.all,
                                  GPS_m,
                                  params = list(alpha = c(0.1,0.2,0.4),
                                                beta=0.2,
                                                g_sigma = 1,
                                                tune_app = "all"))


  # estimate_cerf_gp returns S3 class
  expect_s3_class(cerf_gp_obj)

})
