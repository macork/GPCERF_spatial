#' @title
#' Estimate the Standard Deviation of the Nugget Term in Full Gaussian Process
#'
#' @description
#' Estimates the standard deviations of the nugget term in full GP by
#' calculating the standard deviations of the residuals.
#'
#' @param data A data.frame of observation data.
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#' @param sigma_obs Covariance matrix between observed covariates.
#' @param inv_sigma_obs Inverse of the covariance matrix between observed covariates.
#'
#' @return
#' A scalar of estimated standard deviation of the nugget term in full GP.
#'
#' @export
#'
#'
estimate_noise_gp <- function(data, sigma_obs, inv_sigma_obs){

  noise <- sd(data$Y - arma_mm(sigma_obs - diag(nrow(sigma_obs)),
                               arma_mm(inv_sigma_obs, data$Y)))

  return(noise)
}
