#' @title
#' Calculate Posterior Means for nnGP Model
#'
#' @description
#' Calculates the posterior mean of a point on the CERF based on the nnGP model.
#' This function also returns the weights assigned to all nearest neighbors when
#' calculating the posterior mean.
#'
#' @param hyperparam A set of hyperparameters in the GP model.
#' @param w  A scaler representing the exposure level for the point of interest
#'  on the CERF.
#' @param GPS_w The GPS for all samples when their exposure levels are set
#'  at \code{w}.
#' @param obs_ord A matrix of two columns. First column is the observed
#' exposure levels of all samples; second is the GPS at the observed exposure
#' levels. The rows are in ascending order for the first column.
#' @param y_obs_ord A vector of observed outcome values. The vector is ordered
#' as \code{obs_ord}.
#' @param n_neighbor The number of nearest neighbors on one side
#' (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest
#' neighbors. The total is \code{2*expand*n_neighbor}.
#' @param block_size Number of samples included in a computation block.
#' Mainly used to balance the speed and memory requirement.
#' Larger \code{block_size} is faster, but requires more memory.
#'
#' @return
#' A two-column matrix. The first column is the weights assigned to each
#' nearest neighbor. The second column is the corresponding observed outcome
#' value. The weight in the last row of this matrix is NA and the observed
#' outcome value is the estimated posterior mean of the CERF at point \code{w},
#' which is the weighted sum of all observed outcome values of the neighbors.
#'
#' @export
#'
#' @examples
#'
#' set.seed(1029)
#' data <- generate_synthetic_data(sample_size = 150, gps_spec = 3)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov_mt = as.matrix(data[,-(1:2)]),
#'                    w_all = as.matrix(data$treat))
#'
#' # Hyperparameter
#' hyperparam <- c(0.1, 0.2, 1)
#' n_neighbor <- 10
#' expand <- 1
#' block_size <- 10000
#'
#' # Exposure level
#' wi <- 0.4
#'
#' # Estimate GPS for the exposure level
#' GPS_w <- dnorm(wi,
#'                mean = GPS_m$e_gps_pred,
#'                sd = GPS_m$e_gps_std, log = TRUE)
#'
#' # Order data for easy selection
#' coord_obs = cbind(data$treat, GPS_m$GPS)
#' y_use <- data$Y
#'
#' obs_ord <- coord_obs[order(coord_obs[,1]),]
#' y_use_ord <- y_use[order(coord_obs[,1])]
#'
#' val <- compute_posterior_m_nn(hyperparam = hyperparam,
#'                               w = wi,
#'                               GPS_w = GPS_w,
#'                               obs_ord = obs_ord,
#'                               y_obs_ord = y_use_ord,
#'                               n_neighbor = n_neighbor,
#'                               expand = expand,
#'                               block_size = block_size)
#'
compute_posterior_m_nn <- function(hyperparam,
                                   w,
                                   GPS_w,
                                   obs_ord,
                                   y_obs_ord,
                                   kernel_fn = function(x) exp(-x^2),
                                   n_neighbor = 10,
                                   expand = 5,
                                   block_size = 1e4){


  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  n <- base::length(GPS_w)

  # Compute number of blocks
  n_block <- base::ceiling(n/block_size)


  if(w >= obs_ord[nrow(obs_ord),1]){
    idx_select <- seq( nrow(obs_ord) - expand*n_neighbor + 1, nrow(obs_ord), 1)
  }else{
    idx_anchor <- which.max(obs_ord[,1]>=w)
    idx_start <- max(1, idx_anchor - n_neighbor*expand)
    idx_end <- min(nrow(obs_ord), idx_anchor + n_neighbor*expand)
    if(idx_end == nrow(obs_ord)){
      idx_select <- seq(idx_end - n_neighbor*2*expand + 1, idx_end, 1)
    }else{
      idx_select <- seq(idx_start, idx_start+n_neighbor*2*expand-1, 1)
    }
  }

  used_obs <- t(t(obs_ord[idx_select,])*(1/sqrt(c(alpha, beta))))
  cov_used_inv <- compute_inverse(g_sigma*kernel_fn(as.matrix(dist(used_obs))) + diag(nrow(used_obs)))
  used_y <- y_obs_ord[idx_select]

  w_obs <- t(t(cbind(w, GPS_w))*(1/sqrt(c(alpha, beta))))
  id_all <- split(1:n, ceiling(seq_along(1:n)/n_block))
  all_weights <- sapply(id_all, function(id_ind){
    cov_cross <- g_sigma*kernel_fn(spatstat.geom::crossdist(w_obs[id_ind,1],
                                                            w_obs[id_ind,2],
                                                            used_obs[,1],
                                                            used_obs[,2]))
    #weight
    c(arma_mm(cov_used_inv, Rfast::colsums(cov_cross)))
    # weights_tmp <- cov_cross%*%cov_used_inv
    # w[w<0] <- 0
    # colSums(weights_tmp)
  })
  weights <- Rfast::rowsums(all_weights)/n
  weights[weights<0] <- 0
  weights <- weights/sum(weights)

  est <- c(used_y%*%weights)

  return(cbind(c(idx_select,NA), c(weights, est)))
}
