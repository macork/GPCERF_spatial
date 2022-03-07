#' @title
#' Calculate Posterior Means for nnGP Model
#'
#' @description
#' Calculate the posterior mean of a point on the CERF based on the nnGP model.
#' This function also returns the weights assigned to all nearest neighbours when
#' calculating the posterior mean.
#'
#' @param params Values of hyperparameters in the GP model.
#' @param w.new  The exposure level for the point of interest on the CERF.
#' @param GPS.new The GPS for all samples when their exposure levels are set at \code{w.new}.
#' @param obs.ord A matrix of two columns. First column is the observed exposure levels of all
#' samples; second is the GPS at the observed exposure levels. The rows are in ascending order
#' for the first column.
#' @param y.obs.ord A vector of observed outcome values. The vector is ordered as \code{obs.ord}.
#' @param n.neighbour Number of nearest neighbours on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n.neighbour}.
#' @param block.size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block.size} is faster, but requires more memory.
#'
#' @return
#' A two column matrix. The first column is the weights assigned to each nearest neighbour.
#' The second column is the corresponding observed outcome value. The weight in the last row of
#' this matrix is NA and the observed outcome value is the estimated posterior mean of the CERF
#' at point \code{w.new}, which is the weighted sum of all observed outcome values of the neighbours.
#'
#' @export
#'
#' @examples
get.nn.fast = function(params, w.new, GPS.new, obs.ord, y.obs.ord,
                       n.neighbour = 10, expand = 5, block.size = 1e4){
  # browser()
  n = base::length(GPS.new)
  n.block = base::ceiling(n/block.size)
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
