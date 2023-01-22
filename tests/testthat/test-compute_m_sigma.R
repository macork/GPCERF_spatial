test_that("compute_m_sigma works as expected!", {

   set.seed(1282)
   data <- generate_synthetic_data(sample_size = 250, gps_spec = 3)

   w_all <- seq(0, 20, 0.1)

   #Estimate GPS function
   GPS_m <- train_gps(cov_mt = data[, -(1:2)],
                      w_all = data$treat,
                      sl_lib = c("SL.xgboost"),
                      dnorm_log = FALSE)

   tune_res <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
                               data = data,
                               w = w_all,
                               GPS_m = GPS_m,
                               tuning = FALSE)

   gp_cerf <- tune_res$est

   expect_equal(length(gp_cerf), 201L)
   expect_equal(length(w_all), 201L)
   expect_vector(gp_cerf)
   expect_equal(gp_cerf[10], -26.88945, tolerance = 0.000001)

   tune_res_t <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
                                 data = data,
                                 w = w_all,
                                 GPS_m = GPS_m,
                                 tuning = TRUE)

   expect_equal(length(tune_res_t$cb), 6L)
   expect_equal(length(tune_res_t$est), 201L)
   expect_true(is.na(tune_res_t$est[1]))
})
