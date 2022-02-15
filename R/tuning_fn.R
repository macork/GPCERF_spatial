# tune alpha, beta and gamma in the GP model
tuning.fn = function(param, sim.data, w.all, GPS, e_gps_pred, e_gps_std,
                     kernel.fn = function(x) exp(-x^2)){
  # browser()
  param = unlist(param)
  x.design = model.matrix(~cf1+cf2+cf3+cf4+cf5+cf6-1, data = sim.data)

  obs.use = cbind( sim.data$treat*sqrt(param[1]), GPS*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use))
  inv.Sigma.obs = chol2inv(chol(Sigma.obs))
  # browser()
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
