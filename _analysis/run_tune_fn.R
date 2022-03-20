
librar(xgboost)
# Generate synthetic data
sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)

# Compute true curve
tru.curve <- sapply(seq(0,20,0.1), function(w){
  tru_R(w, sim.data)
})

# Estimate GPS function
e_gps <- xgboost(label=sim.data$treat, data=as.matrix(sim.data[,-(1:2)]), nrounds = 50)
e_gps_pred <- predict(e_gps,as.matrix(sim.data[,-(1:2)]))
e_gps_std <- sd(sim.data$treat-e_gps_pred)
GPS <- dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)

# exposure values
w.all = seq(0,20,0.1)

# tuning + estimating parameters for the GP model
all.params <-  expand.grid(seq(0.01,0.09,0.02), seq(0.01,0.09,0.02), c(1,5,10))
tune.res <- apply(all.params, 1, function(x){
  print(x)
  tuning.fn(param = x, sim.data = sim.data, w.all = w.all, GPS = GPS, e_gps_pred = e_gps_pred,
            e_gps_std = e_gps_std)
})
opt.idx = order(sapply(tune.res, function(x){ mean(x$cb) }))[1]
gp.cerf = tune.res[[opt.idx]]$est


