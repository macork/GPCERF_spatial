
#' @title
#' Estimate the Conditional Exposure Response Function using Gaussian Process
#'
#' @description
#' Estimates the conditional exposure response function (cerf) using Gaussian
#' Process (gp). The function tune the best match (the lowest covariate balance)
#' for the provided set of hyperparameters.
#'
#' @param data A data.table of observation data.
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
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
#' alpha, beta, and g_sigma can be a vector of parameters.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#' @param kernel_fn A kernel function. A default value is a Gaussian Kernel.
#'
#' @return
#' A cerf_gp object that includes the following values:
#'  - TBD
#'
#' @export
#'
#' @examples
#'
#' set.seed(129)
#' sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
#'
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
#'                    w.all = as.matrix(sim.data$treat))
#'
#' # exposure values
#' w.all = seq(0,20,1)
#'
#' data.table::setDT(sim.data)
#'
#' cerf_gp_obj <- estimate_cerf_gp(sim.data,
#'                                 w.all,
#'                                 GPS_m,
#'                                 params = list(alpha = c(0.1),
#'                                               beta=0.2,
#'                                               g_sigma = 1,
#'                                               tune_app = "all"),
#'                                 nthread = 1)
#'
#'
estimate_cerf_gp <- function(data, w, GPS_m, params, nthread = 1,
                             kernel_fn = function(x) exp(-x^2)){


  # Log system info
  log_system_info()

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

  check_params(c("alpha","beta","g_sigma","tune_app"), params)

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

  # Compute m, "confidence interval", and covariate balance for provided
  # hyperparameters. -----------------------------------------------------------
  lfp <- get_options("logger_file_path")

  # make a cluster
  cl <- parallel::makeCluster(nthread, type="PSOCK",
                              outfile= lfp)

  # export variables and functions to cluster cores
  parallel::clusterExport(cl=cl,
                          varlist = c("w", "data", "GPS_m",
                                      "tune_params_subset",
                                      "kernel_fn",
                                      "compute_m_sigma"),
                          envir=environment())

  tune_res <- parallel::parApply(cl, tune_params_subset, 1,
                                 function(x){
                                   compute_m_sigma(hyperparam = x, data = data, w = w,
                                                   GPS_m = GPS_m, kernel_fn = kernel_fn )
  })

  # terminate clusters.
  parallel::stopCluster(cl)

  # Select the combination of hyperparameters that provides the lowest
  # covariate balance ----------------------------------------------------------
  opt_idx <- order(sapply(tune_res, function(x){ mean(x$cb) }))[1]
  gp_cerf <- tune_res[[opt_idx]]$est
  gp_post_sd <- tune_res[[opt_idx]]$pst

  # Build gp_cerf S3 object

  result <- list()
  class(result) <- "cerf_gp"

  result$pst_mean <- gp_cerf
  result$pst_sd <- gp_post_sd
  result$w <- w

  # Add best match to the gp_cerf object

  # Add other useful info form the processing to the gp_cerf object

  # return gp_cerf S3 object
  invisible(result)
}
