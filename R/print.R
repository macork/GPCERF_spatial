#' @title
#' Extend print function for cerf_gp object
#'
#' @param x A cerf_gp object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cerf_gp <- function(x, ...){
  x <- unclass(x)

  cat(" GPCERF Full Gaussian Process exposure rate function object\n")
  cat(" function call: \n")
  cat("      ***       \n")
  print(x$fcall, ...)
  cat("      ***       \n")
  cat(" Look at summary for more details.")
}

#' @title
#' print summary of cerf_gp object
#'
#' @param object A cerf_gp object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data
#' @export
summary.cerf_gp <- function(object, ...){

  cat("GPCERF Full Gaussian Process exposure rate function object summary\n")
  cat_list <- function(input){
    cat(paste("   size: ", length(input),
              ", class: ", class(input),
              ", missing value(s): ", sum(is.na(input)),
              sep = ""))
    if (is.numeric(input)){
      cat(paste("\n   min: ", sprintf("%.3f", min(input, na.rm = TRUE)),
                "\n   max: ", sprintf("%.3f", max(input, na.rm = TRUE)),
                "\n   mean: ", sprintf("%.3f", mean(input, na.rm = TRUE)),
                sep = ""))
    }
  }

  object <- unclass(object)
  object_names <- c("pst_mean", "pst_sd", "w")
  for (item in object_names){
    cat(paste(" ", item, "\n"))
    cat_list(object[[item]])
    cat("\n")
  }
}



#' @title
#' Extend print function for cerf_nngp object
#'
#' @param x A cerf_nngp object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cerf_nngp <- function(x, ...){
  x <- unclass(x)

  cat(" GPCERF Nearest Neighbor Gaussian Process exposure rate function object\n")
  cat(" function call: \n")
  cat("      ***       \n")
  print(x$fcall, ...)
  cat("      ***       \n")
  cat(" Look at summary for more details.")
}



#' @title
#' print summary of cerf_nngp object
#'
#' @param object A cerf_nngp object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data.
#' @export
summary.cerf_nngp <- function(object, ...){

  cat("GPCERF Nearest Neighbore Gaussian Process exposure rate function object summary\n")
  cat_list <- function(input){
    cat(paste("   size: ", length(input),
              ", class: ", class(input),
              ", missing value(s): ", sum(is.na(input)),
              sep = ""))
    if (is.numeric(input)){
      cat(paste("\n   min: ", sprintf("%.3f", min(input, na.rm = TRUE)),
                "\n   max: ", sprintf("%.3f", max(input, na.rm = TRUE)),
                "\n   mean: ", sprintf("%.3f", mean(input, na.rm = TRUE)),
                sep = ""))
    }
  }

  object <- unclass(object)
  object_names <- c("pst_mean", "pst_sd", "w")
  for (item in object_names){
    cat(paste(" ", item, "\n"))
    cat_list(object[[item]])
    cat("\n")
  }
}
