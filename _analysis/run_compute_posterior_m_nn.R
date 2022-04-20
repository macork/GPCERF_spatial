set.seed(3089)
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

# Estimate GPS function
GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                   w.all = as.matrix(data$treat))

# Hyperparameter
hyperparam <- c(0.1, 0.2, 1)
n_neighbor <- 10
expand <- 1
block_size <- 10000

# compute posterior mean and standard deviation for vector of w.
w <- seq(0,20,2)
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])

hyperparam_grid <- expand.grid(seq(0.5,2.5,1),
                               seq(0.4,0.8,0.2),
                               seq(0.5,1.5,1))

optimal_cb <- find_optimal_nn(w_obs = data$treat,
                              w = w,
                              y_obs = data$Y,
                              GPS_m = GPS_m,
                              design_mt = design_mt,
                              hyperparams = hyperparam_grid,
                              n_neighbor = 50, expand = 2, block_size = 2e3)

opt_idx_nn <- order(colMeans(abs(optimal_cb)))[1]
nn_opt_param <- unlist(hyperparam_grid[opt_idx_nn,])

# # Exposure level
# wi <- 0.4
#
# # Estimate GPS for the exposure level
# GPS_w = dnorm(wi,
#               mean = GPS_m$e_gps_pred,
#               sd = GPS_m$e_gps_std, log = T)

# Order data for easy selection
# coord_obs = cbind(data$treat, GPS_m$GPS)
# y_use <- data$Y
#
# obs_ord <- coord_obs[order(coord_obs[,1]),]
# y_use_ord <- y_use[order(coord_obs[,1])]

# # compute noise
# noise <- estimate_noise_nn(hyperparam = hyperparam,
#                            w_obs = data$treat,
#                            GPS_obs = GPS_m$GPS,
#                            y_obs = y_use_ord,
#                            n_neighbor = n_neighbor)

# val <- estimate_mean_sd_nn(hyperparam = hyperparam,
#                                 sigma2 = noise,
#                                 w_obs=data$treat,
#                                 w=w,
#                                 y_obs=y_use_ord,
#                                 GPS_m = GPS_m,
#                                 n_neighbor = 50,
#                                 expand = 1,
#

# compute posterior standard deviation
# pst_sd <- compute_posterior_sd_nn(hyperparam = hyperparam,
#                                   w = wi,
#                                   GPS_w = GPS_w,
#                                   obs_ord = obs_ord,
#                                   sigma2 = noise,
#                                   n_neighbor = 20,
#                                   expand = 1)


# val <- compute_posterior_sd_nn(hyperparam = hyperparam,
#                               w = wi,
#                               GPS_w = GPS_w,
#                               obs_ord = obs_ord,
#                               y_obs_ord = y_use_ord,
#                               n_neighbor = n_neighbor,
#                               expand = expand,
#                               block_size = block_size)
