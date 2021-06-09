library(xgboost)
source("../src/GP_fns.R")
source("../src/sim_fns.R")
source("../src/nnest_fns.R")

args = commandArgs(trailingOnly = T)

sim.data = data_generate(sample_size = as.numeric(args[3]), gps_spec = as.numeric(args[2]))
tru.curve = sapply(seq(0,20,0.1), function(w){
  tru_R(w, sim.data)
})

#estimate GPS functions
e_gps <- xgboost(label=sim.data$treat, data=as.matrix(sim.data[,-(1:2)]), nrounds = 50)
e_gps_pred <- predict(e_gps,as.matrix(sim.data[,-(1:2)]))
e_gps_std <- sd(sim.data$treat-e_gps_pred)
GPS = dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)

# e_gps_pred = (design.mt%*%beta.tru - 0.8)*9+17
# e_gps_std = 5

w.all = seq(0,20,0.1)
# tuning + estimating parameters for the GP model
all.params = expand.grid(seq(0.01,0.09,0.02), seq(0.01,0.09,0.02), c(1,5,10))
tune.res <- apply(all.params, 1, function(x){
  print(x)
  tuning.fn(param = x, sim.data = sim.data, w.all = w.all, GPS = GPS, e_gps_pred = e_gps_pred,
            e_gps_std = e_gps_std)
})
opt.idx = order(sapply(tune.res, function(x){ mean(x$cb) }))[1]
gp.cerf = tune.res[[opt.idx]]$est

## nngp tune+estmate
nn.bc = nn.balance(sim.data$treat, w.all, sim.data$Y, 
           train.GPS.ret = list(GPS = GPS, e_gps_pred = e_gps_pred,e_gps_std_pred = e_gps_std), 
           design.mt = design.mt, all.params = all.params, n.cpu = 1, n.neighbour = 20, expand = 1, block.size = 1)
opt.idx.nn = order(colMeans(abs(nn.bc)))[1]
nn.opt.param = unlist(all.params[opt.idx.nn,])
nn.res = nn.estimate(nn.opt.param, sim.data$treat, w.all, sim.data$Y, 
                      train.GPS.ret = list(GPS = GPS, e_gps_pred = e_gps_pred,e_gps_std_pred = e_gps_std), 
                      n.cpu = 1, n.neighbour = 10, expand = 1)
nn.cerf = sapply(nn.res, function(x) x[nrow(x),2])

# GPS matching
library("GPSmatching")
design.mt = model.matrix(~.-Y-treat-1, data = sim.data)
matched_set = create_matching(sim.data$Y,
                              sim.data$treat,
                              sim.data %>% select(-c(Y, treat)),
                              matching_fun = matching_l1,
                              sl.lib = c("SL.xgboost","SL.earth","SL.gam","SL.ranger"),
                              scale = 0.5,
                              delta_n=1)
matching.cerf = matching_smooth(matched_Y = matched_set$Y,
                      matched_w = matched_set$w,
                      bw.seq = seq(0.2,2,0.2),
                      w.vals = w.all)

# IPTW
library(causaldrf)
iptw.res = iptw_est(Y = Y, treat = treat, treat_formula = treat~cf1+cf2+cf3+cf4+cf5+cf6,
                    numerator_formula = treat~1,
                    data = sim.data,
                    degree = 2,
                    treat_mod = "Normal")
iptw.cerf = iptw.res$param[1] + iptw.res$param[2]*w.all + iptw.res$param[3]*w.all^2

# Kennedy
library(npcausal)
dr.cerf = ctseff(sim.data$Y, sim.data$treat, design.mt, bw.seq = seq(0.2,2,0.2))

saveRDS(list(gp = gp.cerf, nn = nn.cerf, matching = mathcing.cerf, iptw = iptw.cerf, dr = dr.cerf,
             data = sim.data))
