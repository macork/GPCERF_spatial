
t_1 <- proc.time()
set.seed(129)
sim_data <- generate_synthetic_data(sample_size = 200, gps_spec = 1)


# Estimate GPS function
GPS_m <- train_GPS(cov_mt = as.matrix(sim_data[,-(1:2)]),
                   w_all = as.matrix(sim_data$treat))

# exposure values
w_all <- seq(0,20,1)

data.table::setDT(sim_data)

cerf_gp_obj <- estimate_cerf_gp(sim_data,
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
