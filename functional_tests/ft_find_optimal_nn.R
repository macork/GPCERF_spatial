
t_1 <- proc.time()
set.seed(89)
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

 # Estimate GPS function
GPS_m <- train_gps(cov_mt = as.matrix(data[,-(1:2)]),
                   w_all = as.matrix(data$treat),
                   sl_lib = "SL.xgboost",
                   dnorm_log = FALSE)

# Hyperparameter
hyperparam <- c(0.1, 0.2, 1)
n_neighbor <- 10
block_size <- 10000

# compute posterior mean and standard deviation for vector of w.
w <- seq(0,20,2)
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
design_mt <- as.data.frame(design_mt)

hyperparam_grid <- expand.grid(seq(0.5,1.0,1),
                               seq(0.4,0.6,0.2),
                               seq(0.5))

optimal_cb <- GPCERF:::find_optimal_nn(w_obs = data$treat,
                                       w = w,
                                       y_obs = data$Y,
                                       GPS_m = GPS_m,
                                       design_mt = design_mt,
                                       hyperparams = hyperparam_grid,
                                       n_neighbor = 50,
                                       block_size = 2e3,
                                       nthread = 1)
t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))
