library(GPCERF)
set.seed(781)
sim_data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = 12, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
                      w_all = sim_data$treat,
                      sl_lib = c("m_xgboost", "m_ranger"),
                      dnorm_log = TRUE)

# exposure values
q1 <- stats::quantile(sim_data$treat, 0.05)
q2 <- stats::quantile(sim_data$treat, 0.95)

w_all <- seq(q1, q2, 1)

params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all")


t1 <- proc.time()
cerf_gp_obj <- estimate_cerf_gp(sim_data,
                                w_all,
                                GPS_m,
                                params = params_lst,
                                nthread = 12)
t2 <- proc.time()
print(paste("Wall clock time: ", (t2 - t1)[[3]],
            " seconds."))

summary(cerf_gp_obj)
plot(cerf_gp_obj)

# png("readme_gp.png", width = 12, height = 4, units = "in", res = 300)
# plot(cerf_gp_obj)
# dev.off()


