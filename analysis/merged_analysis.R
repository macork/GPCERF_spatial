library(xgboost)
library(snowfall)
library(spatstat)

source("../src/nnest_fns.R")
source("../src/aux.R")
load("/nfs/home/B/bor158/shared_space/ci3_analysis/GPSmatching/prematch_data.RData")
test2 = prematch_data

cov.matrix.no2 = model.matrix( ~.-ZIP - zip - dead - time_count - 
                                 mortality - no2 - mortality_trans,
                               data = test2)
GPS_mod.no2 <-xgboost(data = cov.matrix.no2, label = test2$no2, nrounds=50)
e_gps_pred.no2 <- predict(GPS_mod.no2,cov.matrix.no2)
e_gps_std_pred.no2 <- sd(test2$no2-e_gps_pred.no2)
GPS.no2 <- dnorm(test2$no2, mean = e_gps_pred.no2, sd = e_gps_std_pred.no2)
test2$GPS.no2 = log(GPS.no2)

w.obs.no2 = cbind(test2$no2, test2$GPS.no2)
w.obs.no2.ord = w.obs.no2[order(w.obs.no2[,1]),]
y.obs.no2.ord = test2$mortality_trans[order(w.obs.no2[,1])]
y.un.obs.no2.ord = test2$mortality[order(w.obs.no2[,1])]

sfInit(parallel = T, cpus = 20)
sfExport("get.nn.fast", "e_gps_pred.no2", "e_gps_std_pred.no2", 
         "w.obs.no2.ord", "y.un.obs.no2.ord")
all.res.no2 = sfLapply(seq(0, 80, length.out = 200), function(w){
  print(w)
  GPS.new = log(dnorm(w, mean = e_gps_pred.no2, sd = e_gps_std_pred.no2))
  get.nn.fast(params = c(1,1,1), w.new = w, 
              GPS.new = GPS.new, obs.ord = w.obs.no2.ord, 
              y.obs.ord = y.un.obs.no2.ord, n.neighbour = 50, expand = 6, block.size = 5e3)
})
sfStop()
saveRDS(all.res.no2, "~/cerf/no2_res.rds")
