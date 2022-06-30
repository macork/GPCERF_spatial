set.seed(429)
# generate data
data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
# generate random weights
weights <- runif(nrow(data))
weights <- weights/sum(weights)
# covariate matrix
design_mt <- model.matrix(~.-1, data = data[, 3:ncol(data)])
cb <- GPCERF:::calc_ac(w = data$treat, X = design_mt, weights=weights)
