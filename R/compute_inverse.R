#' @title
#' Compute Matrix Inverse For a Covariate Matrix
#'
#' @description
#' Computes inverse of a covariate matrix using Choleski decomposition.
#'
#' @param mtrx An n*n covariate matrix
#'
#' @return
#' Inverse matrix
#'
#' @export
#'
#' @examples
#'
#' set.seed(934)
#' A <- runif(10)
#' B <- runif(10)
#' C = cbind(A, B)
#' kernel_fn = function(x) exp(-x^2)
#' D = kernel_fn(as.matrix(dist(C)))
#' inv_sigma_obs <- compute_inverse(D)
#'
compute_inverse <- function(mtrx){

  if (!is.matrix(mtrx)){
    stop(paste0("The input mtrx should be a matrix. ",
                "Current format: ", class(mtrx)[1]))
  }

  if (nrow(mtrx) != ncol(mtrx)){
    stop(paste0("The input mtrx should be a square matrix. ",
                "Current dimension: nrow: ",
                nrow(mtrx), ", ncol: ", ncol(mtrx)))
  }


  inv_mtrx <- chol2inv(chol(mtrx))

  return(inv_mtrx)

}
