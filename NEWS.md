## GPCERF (Current Version)

### Fixed
* 

### Changed

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

### Removed
* 

## GPCERF 0.0.1 (2022-03-31)

### Fixed
* 

### Changed

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
