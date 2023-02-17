rm(list = ls())


set.seed(109)
data <- generate_synthetic_data(sample_size = 1000, gps_spec = 1)

data$cf5 <- as.factor(data$cf5)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}


GPCERF::set_logger(logger_level = "TRACE")

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = data[,-(1:2)],
                   w_all = data$treat,
                   sl_lib = c("m_xgboost"),
                   dnorm_log = TRUE)


# exposure values
w_all <- seq(0,20,1)

t_1 <- proc.time()
cerf_nngp_obj <- estimate_cerf_nngp(data,
                                    w_all,
                                    GPS_m,
                                    params = list(alpha = c(0.1, 0.2, 0.3, 0.4),
                                                  beta = 0.2,
                                                  g_sigma = 1,
                                                  tune_app = "all",
                                                  n_neighbor = 50,
                                                  expand = 1,
                                                  block_size = 1000),
                                    formula = ~ . - 1 - Y - treat,
                                    nthread = 12)

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], " s."))

print(cerf_nngp_obj)
summary(cerf_nngp_obj)
#plot(cerf_nngp_obj)

# # profiling
# profvis::profvis({
#   cerf_nngp_obj <- estimate_cerf_nngp(sim_data,
#                                       w_all,
#                                       GPS_m,
#                                       params = list(alpha = c(0.1),
#                                                     beta = 0.2,
#                                                     g_sigma = 1,
#                                                     tune_app = "all",
#                                                     n_neighbor = 20,
#                                                     expand = 1,
#                                                     block_size = 1e4),
#                                       formula = ~ . - 1 - Y - treat,
#                                       nthread = 1)
# })
