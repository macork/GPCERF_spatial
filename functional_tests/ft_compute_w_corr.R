set.seed(124)
data1 <- generate_synthetic_data(sample_size = 200)
data1[,"cf5"] <- as.factor(data1[, "cf5"])

set.seed(987)
weights1 <- runif(nrow(data1))
val1 <- compute_w_corr(w = data1$treat,
                       confounders = data1[, 3:ncol(data1)],
                       weights = weights1)





val2 <- compute_w_corr_2(w = data$treat,
                         covariate = data1[, 3:ncol(data1)],
                         weight = weights1)
