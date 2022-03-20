#' @title
#' Hyperparameter Tuning in Full Guassian Process (GP)
#'
#' @description
#' Calculates the induced covariate balance associated with one hyperparameter
#' configuration in full GP.
#'
#' @param param A vector of values of hyperparameters.
#' @param sim.data A data frame containing all data including outcome, exposure
#' and covariates.
#' @param w.all A vector of exposure levels at which the CERF is estimated.
#' @param GPS A vector of estimated GPS evaluated at the observed exposure levels.
#' @param e_gps_pred A vector of estimated conditional means of the exposure given covariates
#' for all samples.
#' @param e_gps_std A scalar of estimated conditional standard deviations of the exposure given covariates.
#' @param kernel.fn The covariance function of GP.
#'
#' @return
#' A list containing two elements: 1) a vector of absolute weighted correlation of each
#' covariate to the exposure, which is the metric for covariate balance and 2) the estimated
#' CERF at \code{w.all} based on the hyper-parameter values in \code{param}.
#' @export
#'
#' @examples
#'
#' sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)
#'
#' w.all = seq(0,20,0.1)
#'
#' e_gps <- xgboost(label=sim.data$treat, data=as.matrix(sim.data[,-(1:2)]),
#'  nrounds = 50)
#' e_gps_pred <- predict(e_gps,as.matrix(sim.data[,-(1:2)]))
#' e_gps_std <- sd(sim.data$treat-e_gps_pred)
#' GPS <- dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)
#'
#' tune.res <- tuning.fn(param = c(0.09, 0.09, 10), sim.data = sim.data,
#' w.all = w.all, GPS = GPS, e_gps_pred = e_gps_pred,
#' e_gps_std = e_gps_std)
#'
#' gp.cerf <- tune.res$est
#'
tuning.fn = function(param, sim.data, w.all, GPS, e_gps_pred, e_gps_std,
                     kernel.fn = function(x) exp(-x^2)){

  param = unlist(param)
  x.design = model.matrix(~cf1+cf2+cf3+cf4+cf5+cf6-1, data = sim.data)

  obs.use = cbind( sim.data$treat*sqrt(param[1]), GPS*sqrt(param[2]) )
  Sigma.obs = param[3]*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use))
  inv.Sigma.obs = chol2inv(chol(Sigma.obs))


  col.all = sapply(w.all, function(w){

    weights.final = GP.weights.test(w = w, w.obs = sim.data$treat,
                                    obs.use = obs.use, param = param,
                                    inv.Sigma.obs = inv.Sigma.obs,
                                    e_gps_pred = e_gps_pred,
                                    e_gps_std= e_gps_std)
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
