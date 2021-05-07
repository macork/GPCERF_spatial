# Gaussian processes for the estimation of causal exposure-response curves (GP-CERF)
This repository contains the source code for the implementation of GP-CERF and the code used to perform simulation analyses and data application on CERF estimation for the effect of PM2.5 on all-cause mortality. A brief summary of the files in the repo:
* src/GP_fns.R: main functions for GP model where no nearest-neighbour approximation is made.
* src/nnest_fns.R: main functions for GP model with nearest-neighbour approximation.
* src/sim_fns.R: functions used for simulation studies.
* src/utility.cpp: Rcpp files for some utility function used in nnest_fns.R.
* analysis/simulation_analysis.R: Main script to perform simulation studies.
* analysis/merged_analysis.R: Main script to apply GP-CERF on an environmental health dataset to estimate causal effect of PM2.5 on all-cause mortality. No population-stratification is performed.
* analysis/stratified_analysis.R: Main script to apply GP-CERF on an environmental health dataset to estimate causal effect of PM2.5 on all-cause mortality. The dataset is stratified based on age, gender and race. GP-CERF is applied to each stratum separately.
