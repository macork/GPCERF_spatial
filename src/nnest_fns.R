train.GPS = function(cov.mt, w.all){
  require(xgboost)
  GPS_mod <-xgboost(data = cov.mt, label = w.all, nrounds=50)
  e_gps_pred <- predict(GPS_mod,cov.mt)
  e_gps_std_pred <- sd(w.all-e_gps_pred)
  list(GPS = dnorm(w.all, mean = e_gps_pred, sd = e_gps_std_pred),
       e_gps_pred = e_gps_pred, e_gps_std_pred = e_gps_std_pred)
}

get.nn.fast = function(params, w.new, GPS.new, obs.ord, y.obs.ord, 
                       n.neighbour = 10, expand = 5, block.size = 1e4){
  # browser()
  n = length(GPS.new)
  n.block = ceiling(n/block.size)
  #params: length 3, first scale for w, second scale for GPS, 
  #third scale for exp fn
  if(w.new >= obs.ord[nrow(obs.ord),1]){
    idx.all = seq( nrow(obs.ord) - expand*n.neighbour + 1, nrow(obs.ord), 1)
  }else{
    idx.anchor = which.max(obs.ord[,1]>=w.new)
    idx.start = max(1, idx.anchor - n.neighbour*expand)
    idx.end = min(nrow(obs.ord), idx.anchor + n.neighbour*expand)
    if(idx.end == nrow(obs.ord)){
      idx.all = seq(idx.end - n.neighbour*2*expand + 1, idx.end, 1)
    }else{
      idx.all = seq(idx.start, idx.start+n.neighbour*2*expand-1, 1)
    }
  }
  
  obs.use = t(t(obs.ord[idx.all,])*sqrt(params[1:2]))
  cov.use.inv = chol2inv(chol(params[3]*exp(-as.matrix(dist(obs.use))^2) + diag(nrow(obs.use))))
  y.use = y.obs.ord[idx.all]
  
  obs.new = t(t(cbind(w.new, GPS.new))*sqrt(params[1:2]))
  id.all = split(1:n, ceiling(seq_along(1:n)/n.block))
  all.res = sapply(id.all, function(id.ind){
    cov.cross = params[3]*exp(-spatstat.geom::crossdist(obs.new[id.ind,1], obs.new[id.ind,2],
                                                        obs.use[,1], obs.use[,2]))
    #mean
    w = cov.cross%*%cov.use.inv
    y.sum = sum(w%*%y.use)
    w[w<0] = 0
    c(colSums(w), y.sum)
  })
  cbind(c(idx.all,NA), rowSums(all.res)/n)
}

get.nn.sd = function(params, w.new, GPS.new, obs.ord, n.neighbour = 10, expand = 1){
  #params[4] should be sigma^2
  n = length(GPS.new)
  #params: length 3, first scale for w, second scale for GPS, 
  #third scale for exp fn
  if(w.new >= obs.ord[nrow(obs.ord),1]){
    idx.all = seq( nrow(obs.ord) - expand*n.neighbour + 1, nrow(obs.ord), 1)
  }else{
    idx.anchor = which.max(obs.ord[,1]>=w.new)
    idx.start = max(1, idx.anchor - n.neighbour*expand)
    idx.end = min(nrow(obs.ord), idx.anchor + n.neighbour*expand)
    if(idx.end == nrow(obs.ord)){
      idx.all = seq(idx.end - n.neighbour*2*expand + 1, idx.end, 1)
    }else{
      idx.all = seq(idx.start, idx.start+n.neighbour*2*expand-1, 1)
    }
  }
  
  obs.use = t(t(obs.ord[idx.all,])*sqrt(params[1:2]))
  cov.use.inv = chol2inv(chol(params[4]*(params[3]*exp(-as.matrix(dist(obs.use))^2) + 
                                           diag(nrow(obs.use)))))
  obs.new = t(t(cbind(w.new, GPS.new))*sqrt(params[1:2]))
  
  #within variance
  sigma.sq1 = (1+params[3])*params[4]/n
  
  #cross variance
  cross.cov = params[4]*params[3]*exp(-spatstat.geom::crossdist(obs.new[,1],obs.new[,2],
                                                                obs.use[,1],obs.use[,2])^2)
  
  sigma.sq2 = c(calc_cross(cross.cov, cov.use.inv))/n^2
  sqrt(sigma.sq1 - sigma.sq2)
}

calc.ac = function(w, X, weights){
  w.mean = sum(w*weights)
  w.sd = sqrt(sum((w-w.mean)^2*weights))
  w.trans = (w-w.mean)/w.sd
  
  X.mean = colSums(X*weights)
  X.cov = (t(X) - X.mean)%*%diag(weights)%*%t(t(X)-X.mean)
  X.trans = t(t(solve(chol(X.cov)))%*%(t(X)-X.mean))
  
  c(w.trans%*%diag(weights)%*%X.trans)
}

nn.balance = function(w.obs, w.est, y.obs, train.GPS.ret, design.mt, n.cpu = 20,
                      n.neighbour = 50, expand = 2, block.size = 2e3){
  require(snowfall)
  coord.obs = cbind(w.obs, train.GPS.ret$GPS)
  #get rid of unobserved stratified mortality rate
  coord.obs = coord.obs[!is.na(y.obs),]
  y.use = y.obs[!is.na(y.obs)]
  design.use = design.mt[!is.na(y.obs),]
  
  coord.obs.ord = coord.obs[order(coord.obs[,1]),]
  y.use.ord = y.use[order(coord.obs[,1])]
  design.use.ord = design.use[order(coord.obs[,1]),]
  
  all.params = expand.grid(seq(0.5,4.5,1), seq(0.5,4.5,1), seq(0.5,4.5,1))
  sfInit(parallel = T, cpus = n.cpu)
  sfExport("get.nn.fast", "train.GPS.ret", "coord.obs.ord", "y.use.ord", "calc.ac")
  all.cb = apply(all.params, 1, function(params){
    print(params)
    sfExport("params")
    all.res = sfSapply(w.est, function(w){
      print(w)
      GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred, sd = train.GPS.ret$e_gps_std_pred)
      res = get.nn.fast(params = params, w.new = w, GPS.new = GPS.new, obs.ord = coord.obs.ord, 
                  y.obs.ord = y.use.ord, n.neighbour = 50, expand = 2, block.size = 2e3)
      idx = res[-nrow(res),1]
      weights = res[-nrow(res),2]
      weights = weights/sum(weights)
      calc.ac( coord.obs[idx,1], design.use.ord[idx,], weights = weights)
    })
    #covariate specific balance, averaged over w.est
    rowMeans(all.res)
  })
  sfStop()
  all.cb
}

nn.estimate = function(params, w.obs, w.est, y.obs, train.GPS.ret, n.cpu = 20,
                       n.neighbour = 50, expand = 2, block.size = 2e3){
  require(snowfall)
  coord.obs = cbind(w.obs, train.GPS.ret$GPS)
  #get rid of unobserved stratified mortality rate
  coord.obs = coord.obs[!is.na(y.obs),]
  y.use = y.obs[!is.na(y.obs)]
  
  coord.obs.ord = coord.obs[order(coord.obs[,1]),]
  y.use.ord = y.use[order(coord.obs[,1])]
  
  # w.obs = cbind(test2$ozone, test2$GPS.ozone)
  # w.obs.ozone.ord = w.obs.ozone[order(w.obs.ozone[,1]),]
  # y.obs.ozone.ord = test2$mortality_trans[order(w.obs.ozone[,1])]
  # y.un.obs.ozone.ord = test2$mortality[order(w.obs.ozone[,1])]
  
  sfInit(parallel = T, cpus = n.cpu)
  sfExport("get.nn.fast", "train.GPS.ret", "coord.obs.ord", "y.use.ord", "params")
  all.res = sfLapply(w.est, function(w){
    print(w)
    GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred, sd = train.GPS.ret$e_gps_std_pred)
    get.nn.fast(params = params, w.new = w, GPS.new = GPS.new, obs.ord = coord.obs.ord, 
                y.obs.ord = y.use.ord, n.neighbour = 50, expand = 2, block.size = 2e3)
  })
  sfStop()
  all.res
}

nn.sigma.est = function(params, w.obs, GPS.obs, y.obs, n.neighbour, n.core = 20){
  require(snowfall)
  #params: length 3, first scale for w, second scale for GPS, third scale for exp fn
  obs = cbind(w.obs*sqrt(params[1]), GPS.obs*sqrt(params[2]))
  obs.ord = obs[order(w.obs),]
  y.ord = y.obs[order(w.obs)]
  
  sfInit(parallel = T, cpus = n.core)
  sfExport("params","n.neighbour", "w.obs", "obs.ord", "y.ord")
  all.residuals = sfSapply(1:length(w.obs), function(i){
    i.min = max(i-n.neighbour/2,1)
    if(i.min - 1 + n.neighbour >= length(w.obs)){
      idx.use = (length(w.obs)-n.neighbour + 1):(length(w.obs))
    }else{
      idx.use = i.min:(i.min + n.neighbour -1)
    }
    
    dist.all = params[3]*exp(-as.matrix(dist(obs.ord[c(i,idx.use),]))^2) + diag(n.neighbour+1)
    w = dist.all[1,-1]%*%chol2inv(chol(dist.all[-1,-1]))
    c(w%*%y.ord[idx.use]) - y.ord[i]
  })
  sfStop()
  sigma2 = var(all.residuals)
  return(sigma2)
}
