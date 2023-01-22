set_logger(logger_file_path = "functional_tests/GPCERF.log", logger_level = "DEBUG")

t_1 <- proc.time()
set.seed(19)
sim_data <- generate_synthetic_data(sample_size = 120, gps_spec = 3)
# Estimate GPS function
GPS_m <- train_GPS(cov_mt = as.matrix(sim_data[,-(1:2)]),
                  w_all = as.matrix(sim_data$treat))
# exposure values
w_all <- seq(0,20,2)
data.table::setDT(sim_data)
cerf_nngp_obj <- estimate_cerf_nngp(sim_data,
                                    w_all,
                                    GPS_m,
                                    params = list(alpha = c(0.1),
                                                  beta = 0.2,
                                                  g_sigma = 1,
                                                  tune_app = "all",
                                                  n_neighbor = 20,
                                                  expand = 1,
                                                  block_size = 1e4),
                                    formula = ~ . - 1 - Y - treat,
                                    nthread = 12)

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], " s."))

print(cerf_nngp_obj)
summary(cerf_nngp_obj)

# profiling
profvis::profvis({
  cerf_nngp_obj <- estimate_cerf_nngp(sim_data,
                                      w_all,
                                      GPS_m,
                                      params = list(alpha = c(0.1),
                                                    beta = 0.2,
                                                    g_sigma = 1,
                                                    tune_app = "all",
                                                    n_neighbor = 20,
                                                    expand = 1,
                                                    block_size = 1e4),
                                      formula = ~ . - 1 - Y - treat,
                                      nthread = 1)
})
