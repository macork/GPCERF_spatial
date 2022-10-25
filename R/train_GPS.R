#' @title
#' Train A Model for GPS
#'
#' @description
#' Estimates the conditional mean and sd of exposure level as a function of
#' covariates with xgboost algorithm.
#'
#' @param cov_mt A covariate matrix containing all covariates. Each row is a
#' sample and each column is a covariate.
#' @param w_all A vector of observed exposure levels.
#' @param dnorm_log Logical, if TRUE, probabilities p are given as log(p).
#'
#' @return
#' A data.table that includes:
#'   - a vector of estimated GPS at the observed exposure levels;
#'   - a vector of estimated conditional means of exposure levels when the covariates are fixed
#' at the observed values;
#'   - estimated standard deviation of exposure levels
#' @export
#'
#' @examples
#' mydata <- generate_synthetic_data()
#' GPS_m <- train_GPS(mydata[,c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")],
#'                    mydata$treat)
#'

train_GPS <- function(cov_mt, w_all, dnorm_log = FALSE){

  logger::log_info("Started estimating GPS values ... ")
  t_1 <- proc.time()
  # GPS_fit <- CausalGPS::estimate_gps(NA, w_all,
  #                                    cov_mt,
  #                                    pred_model = "sl",
  #                                    gps_model = "parametric",
  #                                    internal_use = T,
  #                                    params = list(xgb_max_depth = c(3,4,5),
  #                                                  xgb_nrounds=c(50,60)),
  #                                    nthread = 1,
  #                                    sl_lib = c("m_xgboost", "m_ranger"))
  #
  # e_gps_pred <- GPS_fit$e_gps_pred
  # e_gps_std <- GPS_fit$e_gps_std_pred
  # GPS <- c(stats::dnorm(w_all, mean = e_gps_pred, sd = e_gps_std, log = dnorm_log))
  # GPS_m <- data.table::data.table(GPS = GPS,
  #                                 e_gps_pred = e_gps_pred,
  #                                 e_gps_std = e_gps_std)
  GPS_SL <- SuperLearner::SuperLearner(Y = w_all, X = cov_mt,
                                      SL.library = c("SL.earth",
                                                     "SL.gam",
                                                     "SL.glm",
                                                     "SL.glm.interaction",
                                                     "SL.mean",
                                                     "SL.ranger"))
  GPS_SL_sd <- sd(w_all - GPS_SL$SL.predict)
  GPS_m <- data.table::data.table(GPS = dnorm(w_all, mean = GPS_SL$SL.predict, sd = GPS_SL_sd, log = dnorm_log),
                         e_gps_pred = GPS_SL$SL.predict,
                         e_gps_std = GPS_SL_sd)


  t_2 <- proc.time()
  logger::log_debug("Wall clock time to estimate GPS values:  ",
                    " {t_2[[3]] - t_1[[3]]} s.")

 return(GPS_m)
}
