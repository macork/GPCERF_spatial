
t_1 <- proc.time()
set.seed(89)
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

hyperparam_grid <- expand.grid(seq(0.5,1.0,1),
                               seq(0.4,0.6,0.2),
                               seq(0.5))

optimal_cb <- find_optimal_nn(w_obs = data$treat,
                              w = w,
                              y_obs = data$Y,
                              GPS_m = GPS_m,
                              design_mt = design_mt,
                              hyperparams = hyperparam_grid,
                              n_neighbor = 50,
                              expand = 2,
                              block_size = 2e3,
                              nthread = 1)
t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))
