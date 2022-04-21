#' @title
#' Find the Optimal Hyper-parameter for the Nearest Neighbor Gaussian Process
#'
#' @description
#' Computes covariate balance for each combincation of provided hyper-parameters
#' and selects the hyper-parameter values that minimizes the covariate balance.
#'
#' @param w_obs A vector of the observed exposure levels.
#' @param w A vector of exposure levels at which CERF will be estimated.
#' @param y_obs A vector of observed outcomes
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param design_mt The covariate matrix of all samples (intercept excluded).
#' @param hyperparams A matrix of candidate values of the hyper-parameters, each row contains a
#' set of values of all hyper-parameters.
#' @param n_neighbor Number of nearest neighbors on one side (see also \code{expand}).
#' @param expand Scaling factor to determine the total number of nearest neighbours. The total is \code{2*expand*n_neighbor}.
#' @param block_size Number of samples included in a computation block. Mainly used to
#' balance the speed and memory requirement. Larger \code{block_size} is faster, but requires more memory.
#'
#' @return
#' Estimated covariate balance scores for the grid of hyper-parameter values considered in \code{hyperparams}.
#' @export
#'
#' @examples
#'
#' set.seed(89)
#' data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
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
#' hyperparam_grid <- expand.grid(seq(0.5,2.5,1),
#'                                seq(0.4,0.6,0.2),
#'                                seq(0.5,1.5,1))
#'
#' optimal_cb <- find_optimal_nn(w_obs = data$treat,
#'                               w = w,
#'                               y_obs = data$Y,
#'                               GPS_m = GPS_m,
#'                               design_mt = design_mt,
#'                               hyperparams = hyperparam_grid,
#'                               n_neighbor = 50, expand = 2, block_size = 2e3)
#'
find_optimal_nn <- function(w_obs, w, y_obs, GPS_m, design_mt,
                      hyperparams = expand.grid(seq(0.5,4.5,1),
                                                seq(0.5,4.5,1),
                                                seq(0.5,4.5,1)),
                      n_neighbor = 50, expand = 2, block_size = 2e3){


  coord.obs = cbind(w_obs, GPS_m$GPS)

  #Remove unobserved outputs
  coord.obs = coord.obs[!is.na(y_obs),]
  y.use = y_obs[!is.na(y_obs)]
  design.use = design_mt[!is.na(y_obs),]

  coord.obs.ord = coord.obs[order(coord.obs[,1]),]
  y.use.ord = y.use[order(coord.obs[,1])]
  design.use.ord = design.use[order(coord.obs[,1]),]


  all.cb = apply(hyperparams, 1, function(hyperparam){
    print(hyperparam)
    all.res = sapply(w, function(wi){
      print(wi)
      # Estimate GPS for requested w.
      GPS_w = dnorm(wi,
                    mean = GPS_m$e_gps_pred,
                    sd = GPS_m$e_gps_std, log = TRUE)

      # Compute posterior mean
      res = compute_posterior_m_nn(hyperparam = hyperparam,
                                   w = wi,
                                   GPS_w = GPS_w,
                                   obs_ord = coord.obs.ord,
                                   y_obs_ord = y.use.ord,
                                   n_neighbor = n_neighbor,
                                   expand = expand,
                                   block_size = block_size)
      idx = res[-nrow(res),1]
      weights = res[-nrow(res),2]
      weights = weights/sum(weights)
      calc_ac( coord.obs[idx,1], design.use.ord[idx,], weights = weights)
    })
    #covariate specific balance, averaged over w
    rowMeans(all.res)
  })

  return(all.cb)

}
