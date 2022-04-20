sim.data = GPCERF::generate_synthetic_data(sample_size = 200)

e_gps <- lm(treat~.-Y, data = sim.data)
e_gps_pred <- e_gps$fitted.values
e_gps_std <- sd(e_gps$residuals)
GPS = dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)

deriv.est = sapply(seq(2.5,17.5,0.1), function(w){
  deriv.nn.fast(w = w, w.obs = sim.data$treat, GPS.obs = GPS, y.obs = sim.data$Y,
              params = c(1,1,1), e_gps_pred = e_gps_pred, e_gps_std = e_gps_std,
              n.neighbour = 20, expand = 1, block.size = 20)
})
