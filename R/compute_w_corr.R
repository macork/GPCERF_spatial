#' @title
#' Compute weighted correlation
#'
#' @description
#' Computes weighted correlation of the observational data based on weights
#' achieved by Gaussian Process.
#'
#' @param w A vector of exposure values for the observed data.
#' @param confounders A data.frame of observational confounders.
#' @param weights A vector of weights for each observation data.
#'
#' @return
#' A vector of covariate balance.
#'
#' @export
#'
#' @examples
#'
#' set.seed(124)
#' mydata <- generate_synthetic_data(sample_size = 200)
#' weights <- runif(nrow(mydata))
#' compute_w_corr(mydata$treat,
#'                mydata[, 3:ncol(mydata)],
#'                weights)
#'
compute_w_corr <- function(w, confounders, weights) {


  logger::log_trace("Computing covariate balance ... ")

  if (!is.data.frame(confounders)) {
    stop(paste0("The confounders should be a data.frame. ",
                "Current format: ", class(confounders)[1]))
  }

  if (!is.vector(w)) {
    stop(paste0("The w param should be a vector. ",
                "Current format: ", class(w)[1]))
  }

  if (nrow(confounders) != length(weights)){
    stop(paste0("Number of data samples in confounder (", nrow(confounders),
                ") and length of ", "weights (", length(weights),
                ") should be equal."))
  }

  if (length(w) != length(weights)){
    stop(paste0("Number of data samples in w (", length(w),
                ") and length of ", "weights (", length(weights),
                ") should be equal."))
  }

  # TODO: model.matrix will create dummy variables for factors.
  # Double-check.
  conf_names <- colnames(confounders)
  frml <- paste("~", paste(conf_names, collapse = "+"), "-1", sep = "")

  # normalize weights here
  weights[weights < 0] <- 0
  if (sum(weights) > 0) {
    weights <- weights / sum(weights)
  }

  x_design <- model.matrix(as.formula(frml), data = confounders)
  w_mean <- sum(w * weights)
  w_sd <- sqrt(sum((w - w_mean) ^ 2 * weights))
  w_stan <- (w - w_mean) / w_sd

  x_mean <- colSums(x_design * weights)
  x_cov <- (t(x_design) - x_mean) %*% diag(weights) %*% t(t(x_design) - x_mean)

  # when x_cov is rank deficient, return NA for all covariate balance
  x_stan <- tryCatch(t(t(solve(chol(x_cov))) %*% (t(x_design) - x_mean)),
                     error = function(e) NA)
  if (!is.na(x_stan[1])) {
    covariate_balance <- abs(c(t(x_stan) %*% diag(weights) %*% w_stan))
  } else {
    covariate_balance = rep(NA, nrow(x_cov))
  }
  return(covariate_balance)
}
