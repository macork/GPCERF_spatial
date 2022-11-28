#' @title
#' Find the Optimal Hyper-parameter for the Nearest Neighbor Gaussian Process
#'
#' @description
#' Computes covariate balance for each combination of provided hyper-parameters
#' and selects the hyper-parameter values that minimizes the covariate balance.
#'
#' @param w_obs A vector of the observed exposure levels.
#' @param w A vector of exposure levels at which CERF will be estimated.
#' @param y_obs A vector of observed outcomes
#' @param GPS_m A data.frame of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample
#'   (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param design_mt The covariate matrix of all samples (intercept excluded).
#' @param hyperparams A matrix of candidate values of the hyper-parameters,
#' each row contains a set of values of all hyper-parameters.
#' @param n_neighbor The number of nearest neighbors on one side
#' (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest
#' neighbors. The total is \code{2*expand*n_neighbor}.
#' @param block_size The number of samples included in a computation block.
#' Mainly used to balance the speed and memory requirement. Larger
#' \code{block_size} is faster, but requires more memory.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#'
#' @return
#' Estimated covariate balance scores for the grid of hyper-parameter values
#' considered in \code{hyperparams}.
#'
#' @export
#'
#' @examples
#'
#' set.seed(89)
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov_mt = data[,-(1:2)], w_all = data$treat)
#'
#' # Hyperparameter
#' hyperparam <- c(0.1, 0.2, 1)
#' n_neighbor <- 10
#' expand <- 1
#' block_size <- 10000
#'
#' # compute posterior mean and standard deviation for vector of w.
#' w <- seq(0,20,2)
#' design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
#'
#' hyperparam_grid <- expand.grid(seq(0.5,1.0,1),
#'                                seq(0.4,0.6,0.2),
#'                                seq(0.5))
#'
#' optimal_cb <- find_optimal_nn(w_obs = data$treat,
#'                               w = w,
#'                               y_obs = data$Y,
#'                               GPS_m = GPS_m,
#'                               design_mt = design_mt,
#'                               hyperparams = hyperparam_grid,
#'                               n_neighbor = 50, expand = 2, block_size = 2e3,
#'                               nthread = 1)
#'
find_optimal_nn <- function(w_obs, w, y_obs, GPS_m, design_mt,
                      hyperparams = expand.grid(seq(0.5,4.5,1),
                                                seq(0.5,4.5,1),
                                                seq(0.5,4.5,1)),
                      n_neighbor = 50, expand = 2, block_size = 2e3,
                      nthread = 1){

  logger::log_info("Started finding optimal values ... ")
  t_opt_1 <- proc.time()

  coord_obs <- cbind(w_obs, GPS_m$GPS)

  #Remove unobserved outputs
  coord_obs <- coord_obs[!is.na(y_obs),]
  y_use <- y_obs[!is.na(y_obs)]
  design_use <- design_mt[!is.na(y_obs),]

  coord_obs_ord <- coord_obs[order(coord_obs[,1]),]
  y_use_ord <- y_use[order(coord_obs[,1])]
  design_use_ord <- design_use[order(coord_obs[,1]),]

  lfp <- get_options("logger_file_path")

  # make a cluster
  t_cl_1 <- proc.time()
  cl <- parallel::makeCluster(nthread, type="PSOCK",
                              outfile= lfp)

  # export variables and functions to cluster cores
  parallel::clusterExport(cl=cl,
                          varlist = c("w", "GPS_m",
                                      "coord_obs_ord", "y_use_ord",
                                      "n_neighbor", "expand", "block_size",
                                      "compute_posterior_m_nn", "calc_ac"),
                          envir=environment())

  t_cl_2 <- proc.time()

  logger::log_debug("Time to setup cluster with {nthread} core(s):",
                   "{t_cl_2[[3]] - t_cl_1[[3]]} s.")

  all_cb <- apply(hyperparams, 1, function(hyperparam){


    # export apply related parameters.
    parallel::clusterExport(cl=cl,
                            varlist = c("hyperparam"),
                            envir=environment())

    all_res_list <- parallel::parLapply(cl,
                                        w,
                                        function(wi){
      # Estimate GPS for requested w.
      GPS_w <- dnorm(wi,
                     mean = GPS_m$e_gps_pred,
                     sd = GPS_m$e_gps_std, log = TRUE)

      # Compute posterior mean
      res <- compute_posterior_m_nn(hyperparam = hyperparam,
                                    w = wi,
                                    GPS_w = GPS_w,
                                    obs_ord = coord_obs_ord,
                                    y_obs_ord = y_use_ord,
                                    n_neighbor = n_neighbor,
                                    expand = expand,
                                    block_size = block_size)
      idx <- res[-nrow(res),1]
      weights <- res[-nrow(res),2]
      weights <- weights/sum(weights)
      calc_ac( coord_obs[idx,1], design_use_ord[idx,], weights = weights)
    })

    all_res <- do.call(cbind, all_res_list)

    #covariate specific balance, averaged over w
    rowMeans(all_res, na.rm = T)
  })

  # terminate clusters.
  parallel::stopCluster(cl)

  t_opt_2 <- proc.time()
  logger::log_info("Done with finding optimal value.",
                   "(Wall clock time: {t_opt_2[[3]] - t_opt_1[[3]]} s.} ... ")


  return(all_cb)
}
