#' @title
#' Title
#'
#' @description
#'
#' @param params
#' @param w.new
#' @param GPS.new
#' @param obs.ord
#' @param n.neighbour
#' @param expand
#'
#' @return
#' @export
#'
#' @examples
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
