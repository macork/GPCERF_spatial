#' @title
#' Hyper-parameter Tuning of nnGP
#'
#' @description
#' Select the optimal hyper-parameter values on a pre-defined grid by minimizing the covariate
#' balance of the nnGP.
#'
#' @param w.obs A vector of the observed exposure levels.
#' @param w.est A vector of exposure levels at which CERF will be estimated.
#' @param y.obs A vector of observed outcomes
#' @param train.GPS.ret A list containing three items, 1) estimated GPS at observed exposure
#' levels; 2) estimated conditional means of the exposure level at the observed covariate
#' values for all samples; 3) estimated conditional standard deviation of the exposure level
#' given all covariates.
#' @param design.mt The covariate matrix of all samples (intercept excluded).
#' @param all.params A matrix of candidate values of the hyper-parameters, each row contains a
#' set of values of all hyper-parameters.
#' @param n.cpu Number of cores to use in the tuning process.
#' @param n.neighbour Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n.neighbour}.
#' @param block.size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block.size} is faster, but requires more memory.
#'
#' @return
#' Estimated covariate balance scores for the grid of hyper-parameter values considered in \code{all.params}.
#' @export
#'
#' @examples
nn.balance = function(w.obs, w.est, y.obs, train.GPS.ret, design.mt,
                      all.params = expand.grid(seq(0.5,4.5,1), seq(0.5,4.5,1), seq(0.5,4.5,1)),
                      n.cpu = 20,  n.neighbour = 50, expand = 2, block.size = 2e3){
  require(snowfall)
  coord.obs = cbind(w.obs, train.GPS.ret$GPS)
  #get rid of unobserved stratified mortality rate
  coord.obs = coord.obs[!is.na(y.obs),]
  y.use = y.obs[!is.na(y.obs)]
  design.use = design.mt[!is.na(y.obs),]

  coord.obs.ord = coord.obs[order(coord.obs[,1]),]
  y.use.ord = y.use[order(coord.obs[,1])]
  design.use.ord = design.use[order(coord.obs[,1]),]

  sfInit(parallel = T, cpus = n.cpu)
  sfExport("get.nn.fast", "train.GPS.ret", "coord.obs.ord", "y.use.ord", "calc.ac")
  all.cb = apply(all.params, 1, function(params){
    print(params)
    sfExport("params")
    all.res = sfSapply(w.est, function(w){
      print(w)
      GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred, sd = train.GPS.ret$e_gps_std_pred, log = T)
      res = get.nn.fast(params = params, w.new = w, GPS.new = GPS.new, obs.ord = coord.obs.ord,
                        y.obs.ord = y.use.ord, n.neighbour = n.neighbour, expand = expand, block.size = block.size)
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
