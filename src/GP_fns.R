library(spatstat)

GP.weights.calc = function(w, w.obs, obs.use, param, inv.Sigma.obs, e_gps_pred, e_gps_std,
                           kernel.fn = function(x) exp(-x^2)){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std)
  
  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                                            obs.use[,1], obs.use[,2]))
  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  weights.all = Sigma.cross%*%inv.Sigma.obs
  weights.all
}

GP.weights.test = function(w, w.obs, obs.use, param, inv.Sigma.obs, e_gps_pred, e_gps_std, 
                           kernel.fn = function(x) exp(-x^2)){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std)
  
  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                                            obs.use[,1], obs.use[,2]))
  # each row is the weights for all subject for estimate of Y_i(w)
  # each column is the weight of an observed sample (w_i, c_i)
  # system.time(eigenMapMatMult(Sigma.cross, inv.Sigma.obs))
  c((rep(1/length(w.obs),length(w.obs))%*%Sigma.cross)%*%inv.Sigma.obs)
}

# tune alpha, beta and gamma in the GP model
tuning.fn = function(param, sim.data, w.all, GPS, e_gps_pred, e_gps_std, 
                     kernel.fn = function(x) exp(-x^2)){
  param = unlist(param)
  x.design = model.matrix(~cf1+cf2+cf3+cf4+cf5+cf6-1, data = sim.data)
  
  obs.use = cbind( sim.data$treat*sqrt(param[1]), GPS*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use))
  inv.Sigma.obs = chol2inv(chol(Sigma.obs))
  
  col.all = sapply(w.all, function(w){
    
    weights.final = GP.weights.test(w = w, w.obs = sim.data$treat, obs.use = obs.use, param = param,
                    inv.Sigma.obs = inv.Sigma.obs,
                    e_gps_pred = e_gps_pred, e_gps_std= e_gps_std)
    weights.final[weights.final<0] = 0
    weights.final = weights.final/sum(weights.final)
    
    est = sim.data$Y%*%weights.final
    
    # weighted correlation
    # this computes rho_r(w) for each covariate r
    w.mean = sum(sim.data$treat*weights.final)
    w.sd = sqrt(sum((sim.data$treat - w.mean)^2*weights.final))
    w.stan = (sim.data$treat - w.mean)/w.sd
    
    x.mean = colMeans(x.design*weights.final)
    x.cov = (t(x.design) - x.mean)%*%diag(weights.final)%*%t(t(x.design) - x.mean)
    x.stan = t(t(solve(chol(x.cov)))%*%(t(x.design) - x.mean))
    c( abs(c(t(x.stan)%*%diag(weights.final)%*%w.stan)), est )
  })
  # this is vector of average rho_r(w) over the range of w for every r
  list(cb = rowMeans(col.all[1:6,]), est = col.all[7,])
}

# estimate sigma
noise.est = function(param, sim.data, GPS){
  # browser()
  n.sample = nrow(sim.data)
  w.all = rep(sim.data$treat, 2)
  GPS.all = rep(GPS, 2)
  Sigma.all = param[3]*exp(-as.matrix(dist(cbind(w.all*sqrt(param[1]), GPS.all*sqrt(param[2]))))) + 
    diag(2*n.sample)
  
  sd(sim.data$Y - Sigma.all[1:n.sample,-(1:n.sample)]%*%chol2inv(chol(Sigma.all[-(1:n.sample),-(1:n.sample)]))%*%sim.data$Y)
}

# weights for derivatives
# we switch back to dnorm to calc GPS
GP.deriv.weights.calc = function(w, w.obs, GPS.obs, param, e_gps_pred, e_gps_std, 
                                 kernel.fn = function(x) exp(-x^2),
                                 kernel.deriv.fn = function(x, mu, sigma) (x-mu)/sigma^2){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS.new = dnorm(w, mean = e_gps_pred, sd = e_gps_std)
  
  obs.use = cbind( w.obs*sqrt(param[1]), GPS.obs*sqrt(param[2]) )
  obs.new = cbind( w*sqrt(param[1]), GPS.new*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use))
  Sigma.cross = param[3]*kernel.fn(spatstat.geom::crossdist(obs.new[,1], obs.new[,2],
                                                            obs.use[,1], obs.use[,2]))*
    param[1]*outer(w, w.obs, "-") - param[2]*outer(GPS.new, GPS.obs, "-")*
    GPS.new*kernel.deriv.fn(w, e_gps_pred, e_gps_std)
  weights.all = Sigma.cross%*%chol2inv(chol(Sigma.obs))
  weights.all
}
