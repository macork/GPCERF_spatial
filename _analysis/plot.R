sim.data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)

# Estimate GPS function
# In the future, CausalGPS gps estimation will be used.
GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                   w.all = as.matrix(sim.data$treat))

# exposure values
w.all = seq(0,20,0.5)

data.table::setDT(sim.data)
cerf_gp_obj <- estimate_cerf_gp(sim.data,
                                w.all,
                                GPS_m,
                                params = list(alpha = c(0.1,0.2,0.4),
                                              beta=0.2,
                                              g_sigma = 1,
                                              tune_app = "all"))

plt_data = data.frame(w = cerf_gp_obj$w, mean = cerf_gp_obj$pst_mean, sd = cerf_gp_obj$pst_sd)
g <- ggplot2::ggplot(plt_data, aes(x = w, y = mean, ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) +
  ggplot2::geom_abline() + ggplot2::geom_ribbon(fill = "red", alpha = 0.25) +
  ggplot2::theme_bw() + ggplot2::ggtitle("Estimated CERF with credible band") + ggplot2::xlab("Exposure level") +
  ggplot2::ylab("Population average counterfactual outcome")


plot(g)

plt_data = data.frame(w = cerf_gp_obj$w, mean = cerf_gp_obj$pst_mean, sd = cerf_gp_obj$pst_sd)

g <- ggplot2::ggplot(plt_data) +
     ggplot2::geom_ribbon(aes(w, y = mean,
                                 ymin = mean - 1.96*sd,
                                 ymax = mean + 1.96*sd),
                          fill = "blue", alpha = 0.25) +
     ggplot2::geom_line(aes(w,mean), color="blue", size = 1) +
     ggplot2::theme_bw() +
     ggplot2::ggtitle("Estimated CERF with credible band (1.96sd)") +
     ggplot2::xlab("Exposure level") +
     ggplot2::ylab("Population average counterfactual outcome")

aes(x = w, y = mean, ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) +
  ggplot2::geom_abline() + ggplot2::geom_ribbon(fill = "red", alpha = 0.25) +
  ggplot2::theme_bw() + ggplot2::ggtitle("Estimated CERF with credible band") + ggplot2::xlab("Exposure level") +
  ggplot2::ylab("Population average counterfactual outcome")


plot(g)
