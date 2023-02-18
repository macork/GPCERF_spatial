library(GPCERF)
set.seed(781)
sim_data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)
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

cerf_gp_obj <- estimate_cerf_gp(sim_data,
                                w_all,
                                GPS_m,
                                params = list(alpha = c(1, 10, 100),
                                              beta = c(0.01, 0.1, 1, 10, 100),
                                              g_sigma = c(0.001, 0.01, 0.1, 1),
                                              tune_app = "all"),
                                nthread = 12)

summary(cerf_gp_obj)
plot(cerf_gp_obj)

png("readme_gp.png", width = 12, height = 4, units = "in", res = 300)
plot(cerf_gp_obj)
dev.off()


