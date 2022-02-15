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
