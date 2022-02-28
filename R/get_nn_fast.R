#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param params param's description
#' @param w.new  param's description
#' @param GPS.new param's description
#' @param obs.ord param's description
#' @param y.obs.ord param's description
#' @param n.neighbour param's description
#' @param expand param's description
#' @param block.size param's description
#'
#' @return
#' @export
#'
#' @examples
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
  all.weights = sapply(id.all, function(id.ind){
    cov.cross = params[3]*exp(-spatstat.geom::crossdist(obs.new[id.ind,1], obs.new[id.ind,2],
                                                        obs.use[,1], obs.use[,2]))
    #mean
    w = cov.cross%*%cov.use.inv
    w[w<0] = 0
    colSums(w)
  })
  weights = rowSums(all.weights)/n
  weights = weights/sum(weights)

  est = c(y.use%*%weights)

  cbind(c(idx.all,NA), c(weights, est))
}
