#' @title
#' Title
#'
#' @description
#' Description
#'
#' @param cov.mt param's description
#' @param w.all param's description
#'
#' @return
#' @export
#'
#' @examples
train.GPS = function(cov.mt, w.all){
  require(xgboost)
  GPS_mod <-xgboost(data = cov.mt, label = w.all, nrounds=50)
  e_gps_pred <- predict(GPS_mod,cov.mt)
  e_gps_std_pred <- sd(w.all-e_gps_pred)
  list(GPS = dnorm(w.all, mean = e_gps_pred, sd = e_gps_std_pred),
       e_gps_pred = e_gps_pred, e_gps_std_pred = e_gps_std_pred)
}
