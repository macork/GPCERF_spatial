#' @title
#' A helper function for cerf_gp object
#'
#' @description
#' A helper function to plot cerf_gp object using ggplot2 package.
#'
#' @param object A cerf_gp object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#' @export
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cerf_gp <- function(object, ...){


  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i,unlist(dot_args[i],use.names = FALSE))
  }

  # extract data
  tmp_data <- data.frame(w_vals = object$w,
                         mean_vals = object$pst_mean,
                         sd_vals = object$pst_sd)

  g <- ggplot2::ggplot(tmp_data) +
       ggplot2::geom_ribbon(ggplot2::aes(.data$w_vals,
                                y = .data$mean_vals,
                                ymin = .data$mean_vals - 1.96*.data$sd_vals,
                                ymax = .data$mean_vals + 1.96*.data$sd_vals),
                                fill = "blue", alpha = 0.25) +
        ggplot2::geom_line(ggplot2::aes(.data$w_vals,.data$mean_vals),
                           color="blue", size = 1) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("Estimated CERF (gp) with credible band (1.96sd)") +
        ggplot2::xlab("Exposure level") +
        ggplot2::ylab("Population average counterfactual outcome")

  return(g)
}


#' @title
#' Extend generic plot functions for cerf_gp class
#'
#' @description
#' A wrapper function to extend generic plot functions for cerf_gp class.
#'
#' @param x  A cerf_gp object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cerf_gp <- function(x, ...){
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}

#' @title
#' A helper function for cerf_nngp object
#'
#' @description
#' A helper function to plot cerf_nngp object using ggplot2 package.
#'
#' @param object A cerf_nngp object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#' @export
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cerf_nngp <- function(object, ...){

  gg_labs <- NULL
  gg_title <- "Exposure Rate Function"

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i,unlist(dot_args[i],use.names = FALSE))
  }

  # extract data
  tmp_data <- data.frame(w_vals = object$w,
                         mean_vals = object$pst_mean,
                         sd_vals = object$pst_sd)

  g <- ggplot2::ggplot(tmp_data) +
       ggplot2::geom_ribbon(ggplot2::aes(.data$w_vals,
                                y = .data$mean_vals,
                                ymin = .data$mean_vals - 1.96*.data$sd_vals,
                                ymax = .data$mean_vals + 1.96*.data$sd_vals),
                                fill = "#FC4E07", alpha = 0.25) +
       ggplot2::geom_line(ggplot2::aes(.data$w_vals,.data$mean_vals),
                          color="#FC4E07",
                          size = 1) +
       ggplot2::theme_bw() +
       ggplot2::ggtitle("Estimated CERF (nngp) with credible band (1.96sd)") +
       ggplot2::xlab("Exposure level") +
       ggplot2::ylab("Population average counterfactual outcome")

  return(g)
}


#' @title
#' Extend generic plot functions for cerf_nngp class
#'
#' @description
#' A wrapper function to extend generic plot functions for cerf_nngp class.
#'
#' @param x  A cerf_nngp object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cerf_nngp <- function(x, ...){
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}
