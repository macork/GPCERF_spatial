library(GPCERF)
library(ggplot2)

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


time_df <- data.frame(core = numeric(0), wc = numeric(0))

for (i in c(1, 2, 4, 8)){
t1 <- proc.time()
cerf_gp_obj <- estimate_cerf_gp(sim_data,
                                w_all,
                                GPS_m,
                                params = params_lst,
                                nthread = i)
t2 <- proc.time()
time_df <- rbind(time_df, data.frame(core = i, wc = (t2 - t1)[[3]]))
print(paste("Wall clock time: ", (t2 - t1)[[3]],
            " seconds."))
}

# summary(cerf_gp_obj)
# plot(cerf_gp_obj)

ggplot(time_df) + geom_line(aes(log2(core), log2(wc))) + geom_point(aes(log2(core), log2(wc)))
