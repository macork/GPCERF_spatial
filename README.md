| Resource    |  Github Actions      |  Code Coverage  |
| ----------  | -------------------- | --------------- |
| Platforms   | Windows, macOS, Linux|    codecov      |
| R CMD check | [![R build status](https://github.com/NSAPH-Software/GPCERF/workflows/R-CMD-check/badge.svg)](https://github.com/NSAPH-Software/GPCERF/actions) | [![codecov](https://codecov.io/gh/NSAPH-Software/GPCERF/branch/develop/graph/badge.svg?token=066ISL822N)](https://app.codecov.io/gh/NSAPH-Software/GPCERF) |


# Gaussian processes for the estimation of causal exposure-response curves (GP-CERF)

## Summary
Gaussian Process (GP) and nearest neighbor Gaussian Process (nnGP) approachs for nonparametric modeling. 

## Installation

```r
library("devtools")
install_github("NSAPH-Software/GPCERF", ref="develop")
library("GPCERF")
```

## Usage

### GP

```r
library(GPCERF)
set.seed(781)
sim_data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)

n_core <- 1

m_xgboost <- function(nthread = n_core, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = n_core, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
                      w_all = sim_data$treat,
                      sl_lib = c("m_xgboost", "m_ranger"),
                      dnorm_log = TRUE)

# exposure values
q1 <- stats::quantile(sim_data$treat, 0.05)
q2 <- stats::quantile(sim_data$treat, 0.95)

w_all <- seq(q1, q2, 1)

params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all")

cerf_gp_obj <- estimate_cerf_gp(sim_data,
                                w_all,
                                GPS_m,
                                params = params_lst,
                                nthread = n_core)
summary(cerf_gp_obj)
plot(cerf_gp_obj)
```
```
GPCERF full Gaussian grocess exposure response function object

Optimal hyper parameters(#trial: 300): 
  alpha = 12.9154966501488   beta = 12.9154966501488   g_sigma = 0.1

Optimal covariate balance: 
  cf1 = 0.072 
  cf2 = 0.082 
  cf3 = 0.062 
  cf4 = 0.068 
  cf5 = 0.056 
  cf6 = 0.082

Original covariate balance: 
  cf1 = 0.222 
  cf2 = 0.112 
  cf3 = 0.175 
  cf4 = 0.318 
  cf5 = 0.198 
  cf6 = 0.257
            ----***----       
```

<p>
<img src="man/figures/png/readme_gp.png" width="900">
</p>


### nnGP

```r
set.seed(781)
sim_data <- generate_synthetic_data(sample_size = 5000, gps_spec = 1)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = 12, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

# Estimate GPS function
GPS_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
                      w_all = sim_data$treat,
                      sl_lib = c("m_xgboost", "m_ranger"),
                      dnorm_log = TRUE)

# exposure values
q1 <- stats::quantile(sim_data$treat, 0.05)
q2 <- stats::quantile(sim_data$treat, 0.95)

w_all <- seq(q1, q2, 1)


params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all",
                   n_neighbor = 50,
                   block_size = 1e3)

cerf_nngp_obj <- estimate_cerf_nngp(sim_data,
                                    w_all,
                                    GPS_m,
                                    params = params_lst,
                                    nthread = 12)
summary(cerf_nngp_obj)
plot(cerf_nngp_obj)
```

```
GPCERF nearest neighbore Gaussian process exposure response function object summary

Optimal hyper parameters(#trial: 300): 
  alpha = 0.0278255940220712   beta = 0.215443469003188   g_sigma = 0.1

Optimal covariate balance: 
  cf1 = 0.058 
  cf2 = 0.071 
  cf3 = 0.087 
  cf4 = 0.066 
  cf5 = 0.076 
  cf6 = 0.088

Original covariate balance: 
  cf1 = 0.115 
  cf2 = 0.137 
  cf3 = 0.145 
  cf4 = 0.296 
  cf5 = 0.208 
  cf6 = 0.225
            ----***----                    
```

<p>
<img src="man/figures/png/readme_nngp.png" width="900">
</p>

## References

Ren, B., Wu, X., Braun, D., Pillai, N. and Dominici, F., 2021. Bayesian modeling for exposure response curve via gaussian processes: Causal effects of exposure to air pollution on health outcomes. arXiv preprint arXiv:2105.03454.
