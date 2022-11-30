
A <- runif(10)
B <- runif(10)
C = cbind(A, B)
kernel_fn = function(x) exp(-x^2)
D = kernel_fn(as.matrix(dist(C)))
inv_sigma_obs <- GPCERF:::compute_inverse(D)
print(inv_sigma_obs)
