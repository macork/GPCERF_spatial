
rm(list = ls())
t_1 <- proc.time()
set.seed(129)
data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)

data$cf5 <- as.factor(data$cf5)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}


GPCERF::set_logger(logger_level = "TRACE")

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = data[,-(1:2)],
                      w_all = data$treat,
                      sl_lib = c("m_xgboost"),
                      dnorm_log = FALSE)

# exposure values
w_all <- seq(0,20,1)

cerf_gp_obj <- estimate_cerf_gp(data,
                                w_all,
                                GPS_m,
                                params = list(alpha = c(0.1, 0.2, 0.3),
                                              beta = c(0.2, 0.4, 0.6),
                                              g_sigma = c(0.5, 0.8),
                                              tune_app = "all"),
                                nthread = 12)

t_2 <- proc.time()
print(paste("Wall clock time: ", t_2[[3]] - t_1[[3]], "s."))

print(cerf_gp_obj)
summary(cerf_gp_obj)

plot(cerf_gp_obj)
#
#
# library(ggplot2)
# library(data.table)
#
# object <- cerf_gp_obj
#
# balance <- data.frame(original = object$cb_org,
#                       adjusted = object$cb)
#
# balance$covar_label <- row.names(balance)
#
# # sort data.frame based on original data correlation values
# balance <- balance[order(balance$original), ]
# covar_label <- balance$covar_label
# row.names(balance) <- NULL
#
# m_balance <- reshape(balance,
#                      direction = "long",
#                      varying = 1:2,
#                      v.names = "value",
#                      times = c("original", "adjusted"),
#                      idvar = "covar_label")
#
# rownames(m_balance) <- NULL
# m_balance$covariates <- rep(seq(1, n_cov, 1), 2)
# colnames(m_balance)[colnames(m_balance) == "time"] <- "Data"
#
# default_gg_title <- "Covariate Balance"
# default_gg_labs <- list(x = "Absolute Weighted Correlation", y= "Covariates")
#
# color_var <- c("#1E88E5", "#FFC107")
#
# g <- ggplot2::ggplot(data = m_balance,
#                      ggplot2::aes(x=.data$value,
#                                   y=.data$covariates,
#                                   color=.data$Data)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_path() +
#   ggplot2::scale_y_discrete(limit = factor(1:n_cov),labels = covar_label) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
#   ggplot2::labs(x = default_gg_labs$x, y = default_gg_labs$y) +
#   ggplot2::ggtitle(default_gg_title) +
#   ggplot2::scale_color_manual(values = color_var)
#
#
