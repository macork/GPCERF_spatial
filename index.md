# GPCERF <a href="https://NSAPH-Software.github.io/GPCERF"><img src="man/figures/png/logo_place_holder.png" align="right" width="120" /></a>

<!-- badges: start -->
[![R build status](https://github.com/NSAPH-Software/GPCERF/workflows/R-CMD-check/badge.svg)](https://github.com/NSAPH-Software/GPCERF/actions)
[![codecov](https://codecov.io/gh/fasrc/GPCERF/branch/develop/graph/badge.svg?token=066ISL822N)](https://codecov.io/gh/fasrc/GPCERF)
<!-- badges: end -->
 


# Gaussian processes for the estimation of causal exposure-response curves (GP-CERF)

## Summary
Gaussian Process (GP) approach for nonparametric modeling. 

## Installation

```r
library("devtools")
install_github("NSAPH-Software/GPCERF", ref="develop")
library("GPCERF")
```

## Usage


```r
  sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)

  # Estimate GPS function
  # In the future, CausalGPS gps estimation will be used.
  GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                     w.all = as.matrix(sim.data$treat))

  # exposure values
  w.all = seq(0,20,0.1)

  data.table::setDT(sim.data)

  cerf_gp_obj <- estimate_cerf_gp(sim.data,
                                  w.all,
                                  GPS_m,
                                  params = list(alpha = c(0.1,0.2,0.4),
                                                beta=0.2,
                                                g_sigma = 1,
                                                tune_app = "all"))

```

## References

Ren, B., Wu, X., Braun, D., Pillai, N. and Dominici, F., 2021. Bayesian modeling for exposure response curve via gaussian processes: Causal effects of exposure to air pollution on health outcomes. arXiv preprint arXiv:2105.03454.
