
t_1 <- proc.time()
set.seed(912)
data <- generate_synthetic_data(sample_size = 2000, gps_spec = 3)

w_all = seq(0,20,1)


#Estimate GPS function
GPS_m <- train_gps(cov_mt = data[,-(1:2)],
                   w_all = data$treat,
                   sl_lib = c("SL.xgboost"),
                   dnorm_log = FALSE)

tune_res <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
                            data = data,
                            w = w_all,
                            GPS_m = GPS_m,
                            tuning = FALSE)

gp_cerf <- tune_res$est

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))
