
t_1 <- proc.time()
set.seed(129)
sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)


# Estimate GPS function
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))

# exposure values
w.all = seq(0,20,1)

data.table::setDT(sim.data)

cerf_gp_obj <- estimate_cerf_gp(sim.data,
                                w.all,
                                GPS_m,
                                params = list(alpha = c(0.1),
                                              beta=0.2,
                                              g_sigma = 1,
                                              tune_app = "all"),
                                nthread = 1)

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))
