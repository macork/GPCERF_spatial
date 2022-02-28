#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param w.obs param's description
#' @param w.est param's description
#' @param y.obs param's description
#' @param train.GPS.ret param's description
#' @param design.mt param's description
#' @param all.params param's description
#' @param n.cpu param's description
#' @param n.neighbour param's description
#' @param expand param's description
#' @param block.size param's description
#'
#' @return
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
