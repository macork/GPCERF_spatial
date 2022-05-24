library(GPCERF)
library(ggplot2)
# exposure values
w.all = seq(0,20,0.1)

#first start with small sample size
sim.data <- generate_synthetic_data(sample_size = 100, gps_spec = 1)
# get the true CERF
ture_curve = sapply(w.all, tru_R, sim.data)

# Estimate GPS function
# In the future, CausalGPS gps estimation will be used.
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))

data.table::setDT(sim.data)

# not run
# error during cholesky
system.time( fit_res <- estimate_cerf_gp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                                              params = list(alpha = 10^seq(-2,0,length.out = 5),
                                                            beta = 10^seq(-2,0,length.out = 5),
                                                            g_sigma = c(0.1,1,10),
                                                            tune_app = "all")) )

system.time( fit_res_nn <- estimate_cerf_nngp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                 params = list(alpha = 10^seq(-2,0,length.out = 5),
                               beta = 10^seq(-2,0,length.out = 5),
                               g_sigma = c(0.1,1,10),
                               tune_app = "all",
                               n_neighbor = 10,
                               expand = 2,
                               block_size = 50)) )


# medium sample size
# might still fail sometimes (cholesky)
sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 1)
# get the true CERF
true_curve = sapply(w.all, tru_R, sim.data)

# Estimate GPS function
# In the future, CausalGPS gps estimation will be used.
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))

data.table::setDT(sim.data)

# running time and resulting estimated CERF
# multi-threading is working as expected
system.time( fit_res <- estimate_cerf_gp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                                         params = list(alpha = 10^seq(0,2,length.out = 5),
                                                       beta = 10^seq(0,2,length.out = 5),
                                                       g_sigma = c(0.1,1,10))) )

system.time( fit_res_nn <- estimate_cerf_nngp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                                              params = list(alpha = 10^seq(-2,0,length.out = 5),
                                                            beta = 10^seq(-2,0,length.out = 5),
                                                            g_sigma = c(0.1,1,10),
                                                            tune_app = "all",
                                                            n_neighbor = 10,
                                                            expand = 2,
                                                            block_size = 50)) )

plot(fit_res) + geom_line(data = data.frame(w = w.all, cerf = true_curve), aes(x = w, y = cerf, color = "True CERF"))
plot(fit_res_nn) plot(fit_res) + geom_line(data = data.frame(w = w.all, cerf = true_curve), aes(x = w, y = cerf, color = "True CERF"))

# large sample size
# might still fail sometimes (cholesky)
sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)
# get the true CERF
true_curve = sapply(w.all, tru_R, sim.data)

# Estimate GPS function
# In the future, CausalGPS gps estimation will be used.
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))

data.table::setDT(sim.data)

# running time and resulting estimated CERF
# multi-threading is working as expected
system.time( fit_res <- estimate_cerf_gp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                                         params = list(alpha = 10^seq(0,2,length.out = 5),
                                                       beta = 10^seq(0,2,length.out = 5),
                                                       g_sigma = c(0.1,1,10))) )

system.time( fit_res_nn <- estimate_cerf_nngp(data = sim.data, w = w.all, GPS_m = GPS_m, nthread = 4,
                                              params = list(alpha = 10^seq(-2,0,length.out = 5),
                                                            beta = 10^seq(-2,0,length.out = 5),
                                                            g_sigma = c(0.1,1,10),
                                                            tune_app = "all",
                                                            n_neighbor = 10,
                                                            expand = 2,
                                                            block_size = 50)) )

plot(fit_res) + geom_line(data = data.frame(w = w.all, cerf = true_curve), aes(x = w, y = cerf, color = "True CERF"))
plot(fit_res_nn) + geom_line(data = data.frame(w = w.all, cerf = true_curve), aes(x = w, y = cerf, color = "True CERF"))
