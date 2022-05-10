test_that("compute_m_sigma works as expected!", {

   set.seed(1282)
   data <- generate_synthetic_data(sample_size = 250, gps_spec = 3)

   w_all <- seq(0,20,0.1)

   data.table::setDT(data)

   #Estimate GPS function
   GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                      w.all = as.matrix(data$treat))

   tune_res <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
                               data = data,
                               w = w_all,
                               GPS_m = GPS_m)

   gp_cerf <- tune_res$est

   expect_equal(length(gp_cerf), 201L)
   expect_equal(length(w_all), 201L)
   expect_vector(gp_cerf)
   expect_equal(gp_cerf[10], -7.9525784, tolerance = 0.000001)

})
