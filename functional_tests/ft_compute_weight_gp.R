
set.seed(814)
#Generate synthetic data
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
w_obs <- obs_exposure <- data$treat
# Choose an exposure level to compute CERF
w = 1.8
# Define kernel function
kernel_fn <- function(x) exp(-x^2)
# Estimate GPS function
GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                  w_all = data$treat,
                  sl_lib = c("SL.xgboost"),
                  dnorm_log = FALSE)
GPS <- GPS_m$GPS
# set hyperparameters
hyperparam <- c(0.1, 0.4, 1)
alpha <- hyperparam[1]
beta <- hyperparam[2]
g_sigma <- hyperparam[3]
# Compute scaled observation data and inverse of covariate matrix.
scaled_obs <- cbind(obs_exposure*sqrt(1/alpha), GPS*sqrt(1/beta))
sigma_obs <- g_sigma*kernel_fn(as.matrix(dist(scaled_obs))) + diag(nrow(scaled_obs))
inv_sigma_obs <- compute_inverse(sigma_obs)
weight <- compute_weight_gp(w = w,
                           w_obs = w_obs,
                           scaled_obs = scaled_obs,
                           hyperparam = hyperparam,
                           inv_sigma_obs = inv_sigma_obs,
                           GPS_m = GPS_m,
                           kernel_fn = kernel_fn)
