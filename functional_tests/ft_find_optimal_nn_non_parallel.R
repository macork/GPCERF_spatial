
t_1 <- proc.time()
set.seed(89)
data <- generate_synthetic_data(sample_size = 2000, gps_spec = 3)

 # Estimate GPS function
GPS_m <- train_gps(cov_mt = as.matrix(data[,-(1:2)]),
                   w_all = as.matrix(data$treat),
                   sl_lib = "SL.xgboost",
                   dnorm_log = TRUE)


# compute posterior mean and standard deviation for vector of w.
w <- seq(-10,20,0.4)
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
design_mt <- as.data.frame(design_mt)

hyperparam_grid <- expand.grid(c(0.1, 0.2, 0.3, 0.4),
                               c(0.001, 0.004, 0.01, 0.04),
                               c(0.4, 0.5, 0.6))


w_obs <- data$treat
w <- w
y_obs <- data$Y
GPS_m <- GPS_m
design_mt <- design_mt
hyperparams <- hyperparam_grid
n_neighbor <- 10
expand <- 1
block_size <- 75
nthread <- 1
kernel_fn = function(x) exp(-x ^ 2)

coord_obs <- cbind(w_obs, GPS_m$GPS)

# remove missing values
coord_obs <- coord_obs[!is.na(y_obs), ]
y_use <- y_obs[!is.na(y_obs)]
design_use <- design_mt[!is.na(y_obs), ]

# order data based on exposure level
coord_obs_ord <- coord_obs[order(coord_obs[, 1]), ]
y_use_ord <- y_use[order(coord_obs[, 1])]
design_use_ord <- design_use[order(coord_obs[, 1]), ]


all_res <- apply(hyperparams, 1, function(hyperparam){

  all_res_list <- lapply(w, function(wi){

          # Compute generalized propensity score for all data.
          GPS_w <- dnorm(wi,
                         mean = GPS_m$e_gps_pred,
                         sd = GPS_m$e_gps_std,
                         log = TRUE)

          # estimate posterior mean
          res <- GPCERF:::compute_posterior_m_nn(hyperparam = hyperparam,
                                                 w = wi,
                                                 GPS_w = GPS_w,
                                                 obs_ord = coord_obs_ord,
                                                 y_obs_ord = y_use_ord,
                                                 n_neighbor = n_neighbor,
                                                 kernel_fn = kernel_fn,
                                                 expand = expand,
                                                 block_size = block_size)

          idx <- res[-nrow(res), 1]
          weights <- res[-nrow(res), 2]
          cb <- GPCERF:::compute_w_corr(w = coord_obs_ord[idx, 1],
                                        confounders = design_use_ord[idx, ], weights)
          list(cb = cb, est = res[nrow(res), 2])
  })

  all_cb_tmp <- do.call(cbind, lapply(all_res_list, '[[', 'cb'))
  all_est_tmp <- sapply(all_res_list, '[[', 'est')

  #covariate specific balance, averaged over w
  list( cb = rowMeans(all_cb_tmp, na.rm = TRUE),
        est = all_est_tmp )

})

# Extract the optimum hyperparameters
optimal_cb_res <- all_res
all_cb_res <- sapply(optimal_cb_res, '[[', 'cb')
opt_idx_nn <- order(colMeans(abs(all_cb_res)))[1]
posterior_mean <- optimal_cb_res[[opt_idx_nn]]$est
nn_opt_param <- unlist(hyperparam_grid[opt_idx_nn, ])

