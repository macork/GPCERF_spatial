#' @title
#' Estimation of CERF with nnGP
#'
#' @description
#' Estimate the posterior mean of the CERF at specified exposure levels with nnGP.
#'
#' @param params A vector of hyperparameters for the nnGP.
#' @param w.obs A vector of observed exposure levels.
#' @param w.est A vector of exposure levels at which the CERF is estimated.
#' @param y.obs A vector of observed outcome values.
#' @param train.GPS.ret A list containing three items, 1) estimated GPS at observed exposure
#' levels; 2) estimated conditional means of the exposure level at the observed covariate
#' values for all samples; 3) estimated conditional standard deviation of the exposure level
#' given all covariates.
#' @param n.cpu Number of cores to use in the tuning process.
#' @param n.neighbour Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n.neighbour}.
#' @param block.size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block.size} is faster, but requires more memory.
#'
#' @return
#' A list of matrices, where each matrix is the returned value from \code{get.nn.fast}.
#' @export
#'
#' @examples
nn.estimate = function(params, w.obs, w.est, y.obs, train.GPS.ret, n.cpu = 20,
                       n.neighbour = 50, expand = 2, block.size = 2e3){
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

  snowfall::sfInit(parallel = T, cpus = n.cpu)
  snowfall::sfExport("get.nn.fast", "train.GPS.ret", "coord.obs.ord", "y.use.ord", "params")
  all.res = snowfall::sfLapply(w.est, function(w){
    print(w)
    GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred, sd = train.GPS.ret$e_gps_std_pred, log = T)
    get.nn.fast(params = params, w.new = w, GPS.new = GPS.new, obs.ord = coord.obs.ord,
                y.obs.ord = y.use.ord, n.neighbour = n.neighbour, expand = expand, block.size = block.size)
  })
  snowfall::sfStop()
  all.res
}
