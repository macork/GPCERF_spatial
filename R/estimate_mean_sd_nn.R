#' @title
#' Estimate the CERF with the nnGP Model
#'
#' @description
#' Estimates the posterior mean of the conditional exposure response function
#' at specified exposure levels with nnGP.
#'
#' @param hyperparam A set of hyperparameters for the nnGP.
#' @param sigma2 A scaler representing \code{sigma^2}.
#' @param w_obs A vector of observed exposure levels.
#' @param w A vector of exposure levels at which the CERF is estimated.
#' @param y_obs A vector of observed outcome values.
#' @param GPS_m A data.frame of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample
#'   (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param n_neighbor The number of nearest neighbors on one side
#' (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest
#' neighbors. The total is \code{2*expand*n.neighbour}.
#' @param block_size The number of samples included in a computation block.
#' Mainly used to balance the speed and memory requirement. Larger
#' \code{block.size} is faster, but requires more memory.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#'
#' @return
#' A list of 2 elements, including:
#' - the returned value from \code{compute_posterior_m_nn}
#' - the returned value from \code{compute_posterior_sd_nn}
#'
#' @keywords internal
#'
estimate_mean_sd_nn <- function(hyperparam,
                                sigma2,
                                w_obs,
                                w,
                                y_obs,
                                GPS_m,
                                n_neighbor = 50,
                                expand = 2,
                                block_size = 2e3,
                                nthread = 1){


  t_est_m_sd_1 <- proc.time()
  logger::log_info("Working on estimating mean and sd using nngp approach ...")


  coord_obs <- cbind(w_obs, GPS_m$GPS)

  #Remove missing outputs
  coord_obs <- coord_obs[!is.na(y_obs),]
  y_use <- y_obs[!is.na(y_obs)]

  coord_obs_ord <- coord_obs[order(coord_obs[,1]),]
  y_use_ord <- y_use[order(coord_obs[,1])]


  lfp <- get_options("logger_file_path")

  # make a cluster
  cl <- parallel::makeCluster(nthread, type="PSOCK",
                              outfile= lfp)


  # install the package on all nodes.
  parallel::clusterEvalQ(cl, {library("GPCERF")})

  # export variables and functions to cluster cores
  parallel::clusterExport(cl=cl,
                          varlist = c("w", "GPS_m", "hyperparam",
                                      "coord_obs_ord", "y_use_ord",
                                      "sigma2",
                                      "n_neighbor", "expand", "block_size",
                                      "compute_posterior_m_nn",
                                      "compute_posterior_sd_nn",
                                      "compute_inverse", "calc_cross"),
                          envir=environment())


  all_res_mean <- parallel::parLapply(cl,
                                      w,
                                      function(wi){
    GPS_w <- dnorm(wi,
                  mean = GPS_m$e_gps_pred,
                  sd = GPS_m$e_gps_std, log = TRUE)

    compute_posterior_m_nn(hyperparam = hyperparam,
                           w = wi,
                           GPS_w = GPS_w,
                           obs_ord = coord_obs_ord,
                           y_obs_ord = y_use_ord,
                           n_neighbor = n_neighbor,
                           expand = expand,
                           block_size = block_size)

  })

  # TODO: repeating GPS_w, merge them together.
  all_res_sd <- parallel::parLapply(cl,
                                    w,
                                    function(wi){
    GPS_w <- dnorm(wi,
                   mean = GPS_m$e_gps_pred,
                   sd = GPS_m$e_gps_std, log = TRUE)

    val <- compute_posterior_sd_nn(hyperparam = hyperparam,
                                   w = wi,
                                   GPS_w = GPS_w,
                                   obs_ord = coord_obs_ord,
                                   sigma2 = sigma2,
                                   n_neighbor = n_neighbor,
                                   expand = expand)

    return(val)
  })

  # terminate clusters.
  parallel::stopCluster(cl)

  t_est_m_sd_2 <- proc.time()

  logger::log_info("Done with estimating mean and sd using nngp approach ",
                   "Wall clock time: {t_est_m_sd_2[[3]] - t_est_m_sd_1[[3]]} s.")


  return(list(all_res_mean, unlist(all_res_sd)))
}
