#' @title
#' Estimate the Conditional Exposure Response Function using Nearest Neighbor Gaussian Process
#'
#' @description
#' Estimates the conditional exposure response function (cerf) using
#' the nearest neighbor (nn) Gaussian Process (gp). The function tune the best
#' match (the lowest covariate balance) for the provided set of hyperparameters.
#'
#' @param data A data.table of observation data.
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#'
#' @param w A vector of exposure level to compute CERF.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param params A list of parameters that is required to run the process.
#' These parameters include:
#'   - alpha: A scaling factor for the GPS value.
#'   - beta: A scaling factor for the exposure value.
#'   - g_sigma: A scaling factor for kernel function (gamma/sigma).
#'   - tune_app: A tuning approach. Available approaches:
#'     - all: try all combinations of hyperparameters.
#'   - expand: Scaling factor to determine the total number of nearest neighbors.
#'   The total is \code{2*expand*n.neighbour}.
#'   - n_neighbor: Number of nearest neighbors on one side.
#'   - block_size: Number of samples included in a computation block. Mainly
#'   used to balance the speed and memory requirement. Larger \code{block_size}
#'   is faster, but requires more memory.
#' alpha, beta, and g_sigma can be a vector of parameters.
#' @param kernel_fn A kernel function. A default value is a Gaussian Kernel.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#'
#' @return
#' #' A cerf_nngp object that includes the following values:
#'  - TBD
#'
#' @export
#'
#' @examples
#'
#' set.seed(19)
#' sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov_mt = as.matrix(sim.data[,-(1:2)]),
#'                    w_all = as.matrix(sim.data$treat))
#' # exposure values
#' w.all <- seq(0,20,2)
#' data.table::setDT(sim.data)
#' cerf_nngp_obj <- estimate_cerf_nngp(sim.data,
#'                                     w.all,
#'                                     GPS_m,
#'                                     params = list(alpha = c(0.1),
#'                                                   beta = 0.2,
#'                                                   g_sigma = 1,
#'                                                   tune_app = "all",
#'                                                   n_neighbor = 20,
#'                                                   expand = 1,
#'                                                   block_size = 1e4),
#'                                     nthread = 1)
#'
#'
estimate_cerf_nngp <- function(data, w, GPS_m, params, kernel_fn, nthread = 1){

  # Log system info
  log_system_info()

  # function call
  fcall <- match.call()

  t_nngp_1 <- proc.time()
  logger::log_info("Working on estimating cerf using nngp approach ...")

  # Double-check input parameters ----------------------------------------------
  if (!is.data.table(data)){
    stop(paste0("Data should be a data.table. ",
                "Current format: ", class(data)[1]))
  }

  if (!is.data.table(GPS_m)){
    stop(paste0("The GPS_m should be a data.table. ",
                "Current format: ", class(GPS_m)[1]))
  }

  check_params <- function(my_param, params){
    for (item in my_param){
      if (!is.element(c(item), names(params))){
        stop(paste0("The required parameter, ", item,", is not provided. ",
                    "Current parameters: ", paste(unlist(names(params)),
                                                  collapse = ", ")))
      }
    }
  }

  check_params(c("alpha", "beta", "g_sigma",
                 "tune_app", "n_neighbor",
                 "expand", "block_size"), params)

  # TODO: Check values of parameters, too.

  # Expand the grid of parameters (alpha, beta, g_sigma) -----------------------
  tune_params <-  expand.grid(getElement(params, "alpha"),
                              getElement(params, "beta"),
                              getElement(params, "g_sigma"))

  if (nrow(tune_params) == 0){
    stop(paste("Something went wrong with tuning parameters. ",
               "The expanded grid has not been generated."))
  }

  # Choose subset of tuning parameters based on tuning approach ----------------

  if (getElement(params, "tune_app") == "all"){
    tune_params_subset <- tune_params
  } else if (getElement(params, "tune_app") == "at_random"){
    stop("This approach is not implemented.")
  }

  # Get other parameters -------------------------------------------------------
  n_neighbor <- getElement(params, "n_neighbor")
  expand <- getElement(params, "expand")
  block_size <- getElement(params, "block_size")

  # Search for the best set of parameters --------------------------------------
  design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
  optimal_cb <- find_optimal_nn(w_obs = data[, c(2)][[1]],
                                w = w,
                                y_obs = data[, c(1)][[1]],
                                GPS_m = GPS_m,
                                design_mt = design_mt,
                                hyperparams = tune_params_subset,
                                n_neighbor = n_neighbor,
                                expand = expand,
                                block_size = block_size,
                                nthread = nthread)



  # Extract the optimum hyperparameters
  opt_idx_nn <- order(colMeans(abs(optimal_cb)))[1]
  nn_opt_param <- unlist(tune_params_subset[opt_idx_nn,])

  # Estimate noise -------------------------------------------------------------
  noise_nn <- estimate_noise_nn(hyperparam = nn_opt_param,
                                w_obs = data[, c(2)][[1]],
                                GPS_obs = GPS_m$GPS,
                                y_obs = data[, c(1)][[1]],
                                n_neighbor = n_neighbor)

  # Compute posterior mean and standard deviation ------------------------------
  posterior_vals <- estimate_mean_sd_nn(hyperparam = nn_opt_param,
                                        sigma2 = noise_nn,
                                        w_obs = data[, c(2)][[1]],
                                        w = w,
                                        y_obs = data[, c(1)][[1]],
                                        GPS_m = GPS_m,
                                        n_neighbor = n_neighbor,
                                        expand = expand,
                                        block_size = block_size,
                                        nthread = nthread)

  posterior_mean <- sapply(posterior_vals[[1]], function(x) x[nrow(x),2])
  posterior_sd <- posterior_vals[[2]]


  t_nngp_2 <- proc.time()
  logger::log_info("Done with estimating cerf using nngp approach ",
                   "Wall clock time: {t_nngp_2[[3]] - t_nngp_1[[3]]} s.")


  # Build nngp_cerf S3 object
  result <- list()
  class(result) <- "cerf_nngp"

  result$w <- w
  result$pst_mean <- posterior_mean
  result$pst_sd <- posterior_sd
  result$fcall <- fcall

  invisible(result)
}
