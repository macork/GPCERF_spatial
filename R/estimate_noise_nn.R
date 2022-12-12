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
#' @export
#'
#' @examples
#'
#' set.seed(425)
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
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
#' wi <- 1.2
#'
#' # Estimate GPS for the exposure level
#' GPS_w = dnorm(wi,
#'               mean = GPS_m$e_gps_pred,
#'               sd = GPS_m$e_gps_std, log = TRUE)
#'
#' # Order data for easy selection
#' coord_obs = cbind(data$treat, GPS_m$GPS)
#' y_use <- data$Y
#'
#' obs_ord <- coord_obs[order(coord_obs[,1]),]
#' y_use_ord <- y_use[order(coord_obs[,1])]
#'
#' noise <- estimate_noise_nn(hyperparam = hyperparam,
#'                            w_obs = data$treat,
#'                            GPS_obs = GPS_m$GPS,
#'                            y_obs = y_use_ord,
#'                            n_neighbor = n_neighbor)
#'
estimate_noise_nn <- function(hyperparam,
                              w_obs,
                              GPS_obs,
                              y_obs,
                              n_neighbor,
                              nthread,
                              kern_fn = function(x) exp(-x^2)){


  t_1 <- proc.time()
  logger::log_info("Working on estimating residual error with nngp ...")

  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  obs <- cbind(w_obs*sqrt(1/alpha), GPS_obs*sqrt(1/beta))
  obs_ord <- obs[order(w_obs),]
  y_ord <- y_obs[order(w_obs)]

  lfp <- get_options("logger_file_path")
  cl <- parallel::makeCluster(nthread, type="PSOCK",
                              outfile= lfp)
  parallel::clusterExport(cl=cl,
                          varlist = c("w_obs","obs_ord", "y_ord",
                                      "n_neighbor", "kern_fn", "g_sigma",
                                      "arma_mm"),
                          envir=environment())

  all_residuals <- parallel::parSapply(cl,
                                       1:length(w_obs),
                                       function(i){
    i_min <- max(i-n_neighbor/2,1)
    if(i_min - 1 + n_neighbor >= length(w_obs)){
      idx_use <- (length(w_obs)-n_neighbor + 1):(length(w_obs))
    }else{
      idx_use <- i_min:(i_min + n_neighbor - 1)
    }

    cov_all <- g_sigma*kern_fn(as.matrix(dist(obs_ord[c(i,idx_use),]))) + diag(n_neighbor+1)
    w <- c(arma_mm(compute_inverse(cov_all[-1,-1]), cov_all[1,-1]))
    c(w%*%y_ord[idx_use]) - y_ord[i]
  })
  parallel::stopCluster(cl)

  t_2 <- proc.time()
  logger::log_info("Done with estimating residual error with nngp",
                   "Wall clock time: {t_2[[3]] - t_1[[3]]} s.")

  return( sd(all_residuals) )
}
