#' @title
#' Train Model for GPS
#'
#' @description
#' Estimate the conditional mean and sd of exposure level as a function of covariates with
#' xgboost algorithm.
#'
#' @param cov.mt Covariate matrix containing all covariates. Each row is a sample and each
#' column is a covariate.
#' @param w.all A vector of observed exposure levels.
#'
#' @return
#' A list of three elements: 1) a vector of estimated GPS at the observed exposure levels; 2)
#' a vector of estimated conditional means of exposure levels when the covariates are fixed
#' at the observed values; 3) estimated standard deviation of exposure levels
#' @export
#'
#' @examples
#'
#' mydata <- generate_synthetic_data()
#' gps_list <- train.GPS(as.matrix(mydata[,c("cf1", "cf2", "cf3", "cf4",
#'                                           "cf5", "cf6")]),
#'                       as.matrix(mydata$treat))
#'
#'
train.GPS = function(cov.mt, w.all){
  GPS_mod <- xgboost::xgboost(data = cov.mt, label = w.all, nrounds=50)
  e_gps_pred <- predict(GPS_mod,cov.mt)
  e_gps_std_pred <- sd(w.all-e_gps_pred)
  list(GPS = stats::dnorm(w.all, mean = e_gps_pred, sd = e_gps_std_pred),
       e_gps_pred = e_gps_pred, e_gps_std_pred = e_gps_std_pred)
}
