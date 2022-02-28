#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param params param's description
#' @param w.obs param's description
#' @param w.est param's description
#' @param y.obs param's description
#' @param train.GPS.ret param's description
#' @param n.cpu param's description
#' @param n.neighbour param's description
#' @param expand param's description
#' @param block.size param's description
#'
#' @return
#' @export
#'
#' @examples
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
    GPS.new = dnorm(w, mean = train.GPS.ret$e_gps_pred, sd = train.GPS.ret$e_gps_std_pred, log = T)
    get.nn.fast(params = params, w.new = w, GPS.new = GPS.new, obs.ord = coord.obs.ord,
                y.obs.ord = y.use.ord, n.neighbour = n.neighbour, expand = expand, block.size = block.size)
  })
  sfStop()
  all.res
}
