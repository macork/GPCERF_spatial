library(SuperLearner)

source("../src/GP_fns.R")
source("../src/sim_fns.R")

args = commandArgs(trailingOnly = T)

sim.data = data_generate(sample_size = as.numeric(args[3]), gps_spec = as.numeric(args[2]))
tru.curve = sapply(seq(0,20,0.1), function(w){
  tru_R(w, sim.data)
})

e_gps <- SuperLearner(Y=sim.data$treat, X=sim.data[,-(1:2)], SL.library = c("SL.xgboost","SL.earth","SL.gam","SL.ranger"))
e_gps_pred <- e_gps$SL.predict
e_gps_std <- SuperLearner(Y=abs(sim.data$treat-e_gps_pred), X=sim.data[,-(1:2)], SL.library=c("SL.xgboost","SL.earth","SL.gam","SL.ranger"))
e_gps_std_pred <- e_gps_std$SL.predict
w_resid <- (sim.data$treat-e_gps_pred)/e_gps_std_pred
GPS <- approx(density(w_resid,na.rm = TRUE)$x,density(w_resid,na.rm = TRUE)$y,xout=w_resid,rule=2)$y

w.all = seq(0,20,0.5)
all.params = expand.grid(seq(0.15,0.25,0.01), seq(0.4,1.2,0.1), c(0.5,1,1.5))
system.time(tune.res <- apply(all.params, 1, function(x){
  print(x)
  tuning.fn(param = x, sim.data = sim.data, w.all = w.all, GPS = GPS, e_gps_pred = e_gps_pred,
            e_gps_std_pred = e_gps_std_pred, w_resid = w_resid)
}))

saveRDS(list(data = sim.data, truth = tru.curve, res = tune.res), 
        sprintf("~/cerf/sim_res/scen%s/rep_%s_N_%s.rds", args[2], args[3], args[1]))

#params.opt = unlist(all.params[order(colMeans(tune.res))[1],])
#nugget.est = noise.est(params.opt, sim.data, GPS)
