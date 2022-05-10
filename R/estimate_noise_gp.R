#' @title
#' Estimate the Standard Deviation of the Nugget Term in Full Gaussian Process
#'
#' @description
#' Estimates the standard deviations of the nugget term in full GP by calculating
#' the standard deviations of the residuals.
#'
#' @param hyperparam A vector of hyper-parameter values for the full GP.
#' @param data A data.table of observation data.
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#' @param GPS A vector of estimated GPS at the observed exposure levels.
#'
#' @return
#' A scalar of estimated standard deviation of the nugget term in full GP.
#' @export
#'
#' @examples
#'
#' set.seed(109)
#' data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)
#' data.table::setDT(data)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' hyperparam <- c(0.1, 0.2, 1)
#'
#' noise_est <- estimate_noise_gp(hyperparam, data, GPS_m$GPS)
#'
#'
estimate_noise_gp <- function(hyperparam, data, GPS){


  # Double-check input parameters ----------------------------------------------
  if (!is.data.table(data)){
    stop(paste0("Data should be a data.table. ",
                "Current format: ", class(data)[1]))
  }

  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]

  n_sample <- nrow(data)
  w_all <- rep(data[,2][[1]], 2)
  GPS_all <- rep(GPS, 2)
  Sigma_all <- g_sigma*exp(-as.matrix(dist(cbind(w_all*sqrt(1/alpha),
                                                 GPS_all*sqrt(1/beta))))) +
    diag(2*n_sample)

  noise <- sd(data[,1][[1]] - Sigma_all[1:n_sample,
                        -(1:n_sample)]%*%chol2inv(
                          chol(Sigma_all[-(1:n_sample),
                                         -(1:n_sample)]))%*%data[,1][[1]])

  return(noise)
}
