

set.seed(19)
sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
# Estimate GPS function
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))
# exposure values
w.all <- seq(0,20,0.5)
data.table::setDT(sim.data)
cerf_nngp_obj <- estimate_cerf_nngp(sim.data,
                                    w.all,
                                    GPS_m,
                                    params = list(alpha = c(0.1,0.2),
                                                  beta = 0.2,
                                                  g_sigma = 1,
                                                  tune_app = "all",
                                                  n_neighbor = 20,
                                                  expand = 1,
                                                  block_size = 1e4))


library(ggplot2)

ggplot() + geom_line(aes(cerf_nngp_obj$w, cerf_nngp_obj$pst_mean), color="blue") +
           geom_line(aes(cerf_nngp_obj$w, cerf_nngp_obj$pst_mean - cerf_nngp_obj$pst_sd), color="red") +
           geom_line(aes(cerf_nngp_obj$w, cerf_nngp_obj$pst_mean + cerf_nngp_obj$pst_sd), color="red")


set.seed(429)

# generate data
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

# generate random weights
weights <- runif(nrow(data))
weights <- weights/sum(weights)

# covariate matrix
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])

cb <- calc_ac(w = data$treat, X = design_mt, weights=weights)
