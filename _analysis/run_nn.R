library(xgboost)

source("_analysis/calc_ac.R")
source("_analysis/get_nn_fast.R")

# Step 1: Generate synthetic data
set.seed(591)
sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)

# Step 2: Estimate GPS values
e_gps <- xgboost(label=sim.data$treat, data=as.matrix(sim.data[,-(1:2)]), nrounds = 50)
e_gps_pred <- predict(e_gps,as.matrix(sim.data[,-(1:2)]))
e_gps_std <- sd(sim.data$treat-e_gps_pred)
GPS = dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)

# Step 3: Hyper parameters
w.all = seq(0,20,0.1)
# tuning + estimating parameters for the GP model
all.params = expand.grid(seq(0.01,0.09,0.02), seq(0.01,0.09,0.02), c(1,5,10))

# Step 4: Design matrix

design.mt = model.matrix(~.-Y-treat-1, data = sim.data)

# Step 5: Find the best hyper parameters

## nn_balance function (start) -------------------------------------------------

w.obs <- sim.data$treat
train.GPS.ret <- list(GPS = GPS,
                      e_gps_pred = e_gps_pred,
                      e_gps_std_pred = e_gps_std)
y.obs <-  sim.data$Y
w.est <- w.all
expand <- 1
n.neighbour <- 20
block.size <- 1

coord.obs = cbind(w.obs, train.GPS.ret$GPS)

#get rid of unobserved stratified mortality rate
coord.obs = coord.obs[!is.na(y.obs),]
y.use = y.obs[!is.na(y.obs)]
design.use = design.mt[!is.na(y.obs),]

coord.obs.ord = coord.obs[order(coord.obs[,1]),]
y.use.ord = y.use[order(coord.obs[,1])]
design.use.ord = design.use[order(coord.obs[,1]),]


all.cb = apply(all.params, 1, function(params){
  print(params)
  #sfExport("params")
  all.res = sapply(w.est, function(w){
    print(w)
    # compute GPS for new w with mean and sd of observation
    GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred,
                       sd = train.GPS.ret$e_gps_std_pred,
                       log = T)

    # calculate posterior mean
    res = get.nn.fast(params = params,
                      w.new = w,
                      GPS.new = GPS.new,
                      obs.ord = coord.obs.ord,
                      y.obs.ord = y.use.ord,
                      n.neighbour = n.neighbour,
                      expand = expand,
                      block.size = block.size)

    idx = res[-nrow(res),1]
    weights = res[-nrow(res),2]
    weights = weights/sum(weights)
    calc.ac( coord.obs[idx,1], design.use.ord[idx,], weights = weights)
  })
  #covariate specific balance, averaged over w.est
  rowMeans(all.res)
})

all.cb

## nn_balance function (end) ---------------------------------------------------

# Collect the optimized parameters
nn.bc <- all.cb
opt.idx.nn = order(colMeans(abs(nn.bc)))[1]
nn.opt.param = unlist(all.params[opt.idx.nn,])

## nn_estimate function (start) -------------------------------------------------

params <- nn.opt.param

coord.obs = cbind(w.obs, train.GPS.ret$GPS)
#get rid of unobserved stratified mortality rate
coord.obs = coord.obs[!is.na(y.obs),]
y.use = y.obs[!is.na(y.obs)]

coord.obs.ord = coord.obs[order(coord.obs[,1]),]
y.use.ord = y.use[order(coord.obs[,1])]


all.res = lapply(w.est, function(w){
  print(w)
  GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred,
                  sd = train.GPS.ret$e_gps_std_pred, log = T)
  get.nn.fast(params = params,
              w.new = w,
              GPS.new = GPS.new,
              obs.ord = coord.obs.ord,
              y.obs.ord = y.use.ord,
              n.neighbour = n.neighbour,
              expand = expand,
              block.size = block.size)
})

all.res

## nn_estimate function (end) ---------------------------------------------------

nn.res <- all.res
nn.cerf = sapply(nn.res, function(x) x[nrow(x),2])

plot(nn.cerf)

