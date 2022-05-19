#' @title
#' Calculate True CERF Values
#'
#' @description
#' Derives the true CERF at exposure levels of interest when the data
#' is distributed as described in the simulation study.
#'
#' @param w A exposure level at which the CERF is derived.
#' @param sim.data A data frame from one simulation run.
#'
#' @return
#' The CERF value evaluated at \code{w}.
#' @export
#'
#' @examples
#'
#' sim.data <- generate_synthetic_data()
#' w <- 0.2
#' tr <- tru_R(w, sim.data)
#'
tru_R<-function(w, sim.data){
  design.mt = model.matrix(~cf1+cf2+cf3+cf4+cf5+cf6-1, data = sim.data)
  mean(apply(design.mt, 1, function(x){
    -10 - sum(c(2, 2, 3, -1,2,2)*x) -
      w*(0.1 - 0.1*x[1] + 0.1*x[4] + 0.1*x[5] + 0.1*x[3]^2) + 0.13^2*w^3
  }))
}
