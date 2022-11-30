rm(list = ls())
t_1 <- proc.time()
set.seed(129)
data <- generate_synthetic_data(sample_size = 200, gps_spec = 1)


# Estimate GPS function
GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                   w_all = data$treat,
                   sl_lib = c("SL.xgboost"),
                   dnorm_log = FALSE)

# exposure values
w_all <- seq(0,20,1)

cerf_gp_obj <- estimate_cerf_gp(data,
                                w_all,
                                GPS_m,
                                params = list(alpha = c(0.1),
                                              beta=0.2,
                                              g_sigma = 1,
                                              tune_app = "all"),
                                nthread = 2)

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))

print(cerf_gp_obj)
summary(cerf_gp_obj)

plot(cerf_gp_obj)
