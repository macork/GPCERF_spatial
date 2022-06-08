
t_1 <- proc.time()
set.seed(912)
data <- generate_synthetic_data(sample_size = 250, gps_spec = 3)

w.all = seq(0,20,1)

data.table::setDT(data)

#Estimate GPS function
GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                    w.all = as.matrix(data$treat))

tune_res <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
                            data = data,
                            w = w.all,
                            GPS_m = GPS_m,
                            nthread = 1)

gp.cerf <- tune_res$est

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))
