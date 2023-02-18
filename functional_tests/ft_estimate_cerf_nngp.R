

set.seed(967)
sim_data <- generate_synthetic_data(sample_size = 10000, gps_spec = 1)
sim_data$cf5 <- as.factor(sim_data$cf5)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
                      w_all = sim_data$treat,
                      sl_lib = c("m_xgboost"),
                      dnorm_log = TRUE)

# exposure values
w_all <- seq(min(sim_data$treat)+1, max(sim_data$treat)-1, 1)

cerf_nngp_obj <- estimate_cerf_nngp(sim_data,
                                    w_all,
                                    GPS_m,
                                    params = list(alpha = c(0.01, 0.1, 1,
                                                            2, 3, 4, 8, 16),
                                                  beta = c( 10, 100, 200),
                                                  g_sigma = c(0.0001, 0.001,
                                                              0.01, 0.1),
                                                  tune_app = "all",
                                                  n_neighbor = 100,
                                                  block_size = 1e3),
                                  nthread = 12)




summary(cerf_nngp_obj)
plot(cerf_nngp_obj)

png("readme_nngp.png", width = 8, height = 4, units = "in", res = 300)
plot(cerf_nngp_obj)
dev.off()

