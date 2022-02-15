library(xgboost)
library(snowfall)
library(spatstat)

source("../src/nnest_fns.R")

args = commandArgs(trailingOnly = T)
#args[1] which subpopulation
sub.pop = as.numeric(args[1])

load("/nfs/home/B/bor158/shared_space/ci3_analysis/GPSmatching/prematch_data.RData")
test2 = prematch_data
all.params = expand.grid(sex = 1:2, race = 0:1, age.2 = 0:1)
y.strat = readRDS(sprintf("~/cerf/strat_data/strat_mr_s%s_r%s_a%s.rds", 
                          all.params[sub.pop,1],
                          all.params[sub.pop,2],
                          all.params[sub.pop,3]))

cov.matrix = model.matrix(~.-1 -ZIP - zip - dead - time_count -  mortality - mortality_trans, data = test2)
# 
train.GPS.pm = train.GPS(cov.matrix[,colnames(cov.matrix)!="pm25_ensemble"], test2$pm25_ensemble)
train.GPS.oz = train.GPS(cov.matrix[,colnames(cov.matrix)!="ozone"], test2$ozone)
train.GPS.no2 = train.GPS(cov.matrix[,colnames(cov.matrix)!="no2"], test2$no2)

params.pm = c(0.46,0.32,1.5)
params.oz = c(0.6,0.3,2)
params.no2 = c(0.4,0.3,1)

print("PM2.5")
nn.estimate(params.pm, w.obs = test2$pm25_ensemble, w.est = seq(0,20,length.out = 200), y.obs = y.strat, 
            train.GPS.ret = train.GPS.pm, 
            file.out = sprintf("subpop_res/subpop_%s_pm.rds", sub.pop))

print("Ozone")
nn.estimate(params.oz, w.obs = test2$ozone, w.est = seq(30,50,length.out = 200), y.obs = y.strat, 
            train.GPS.ret = train.GPS.oz, 
            file.out = sprintf("subpop_res/subpop_%s_oz.rds", sub.pop))

print("NO2")
nn.estimate(params.no2, w.obs = test2$no2, w.est = seq(0,60,length.out = 200), y.obs = y.strat, 
            train.GPS.ret = train.GPS.no2, 
            file.out = sprintf("subpop_res/subpop_%s_no2.rds", sub.pop))
