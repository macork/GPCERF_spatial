## GPCERF (development)

### Changed
* find_optimal_nn becomes an internal function. 
* estimate_noise_gp and estimate_noise_nn become an internal function. 
* estimate_mean_sd_nn becomes an internal function.
* compute_weight_gp becomes an internal function.
* compute_w_corr accepts w and confounders separately.  
* compute_posterior_sd_nn becomes an internal function.
* compute_posterior_m_nn becomes an internal function.
* compute_derive_weights_gp becomes an internal function. 
* compute_m_sigma becomes an internal function.
* compute_inverse becomes an internal function. 
* In compute_m_sigma, tuning option does not have a default value. 
* train_gps does not have default values.
* train_gps accepts vector of the SuperLearner package's libraries.
* train_GPS -> train_gps

### Added
* 

## GPCERF 0.1.0 (2022-07-02)

### Changed

* nn_cp_calc -> compute_rl_deriv_nn
* deriv_nn_fast -> compute_deriv_nn
* get_nn_sd -> compute_posterior_sd_nn
* nn_sigma_est -> estimate_noise_nn
* idx.all -> idx_select
* GPS.new -> GPS_w
* w.new -> w
* get.nn.fast -> compute_posterior_m_nn
* w.est -> w 
* nn_balance -> best_nn_cb

### Added

* Package website using pkgdown
* Logger functions
* compute_sd_gp function


## GPCERF 0.0.1 (2022-03-31)

### Fixed
* 

### Changed

* Removed examples from internal functions
* w.obs -> w_obs
* inv.Sigma.obs -> inv_sigma_obs
* obs.use -> scaled_obs
* tune.fn -> compute_m_sigma
* GP.weights.test -> compute_weight_gp
* data.generate -> generate_synthetic_data 


### Added

* estimate_noise function
* estimate_cerf_gp function
* compute_inverse function
* compute_w_corr function

### Removed
* 
