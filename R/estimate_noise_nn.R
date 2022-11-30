#' @title
#' Estimate the Standard Deviation (noise) of the Nugget Term in nnGP
#'
#' @description
#' Estimate the standard deviations of the nugget term (noise) in nnGP by
#' calculating the standard deviations of the residuals.
#'
#' @param hyperparam A vector of hyper-parameter values.
#' @param w_obs A vector of observed exposure levels.
#' @param GPS_obs A vector of estimated GPS evaluated at the observed exposure levels.
#' @param y_obs A vector of observed outcomes.
#' @param n_neighbor Number of nearest neighbors on one side.
#'
#' @return
#' A scalar of estimated standard deviation of the nugget term in nnGP.
#'
#' @keywords internal
#'
estimate_noise_nn <- function(hyperparam,
                              w_obs,
                              GPS_obs,
                              y_obs,
                              n_neighbor){


  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  obs <- cbind(w_obs*sqrt(1/alpha), GPS_obs*sqrt(1/beta))
  obs_ord <- obs[order(w_obs),]
  y_ord <- y_obs[order(w_obs)]

  all_residuals <- sapply(1:length(w_obs), function(i){
    i_min <- max(i-n_neighbor/2,1)
    if(i_min - 1 + n_neighbor >= length(w_obs)){
      idx_use <- (length(w_obs)-n_neighbor + 1):(length(w_obs))
    }else{
      idx_use <- i_min:(i_min + n_neighbor -1)
    }

    dist_all <- g_sigma*exp(-as.matrix(dist(obs_ord[c(i,idx_use),]))^2) + diag(n_neighbor+1)
    w <- dist_all[1,-1]%*%chol2inv(chol(dist_all[-1,-1]))
    c(w%*%y_ord[idx_use]) - y_ord[i]
  })

  return( sd(all_residuals) )
}
