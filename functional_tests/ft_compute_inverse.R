
kernel_fn = function(x) exp(-x^2)

df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c('size', 'wc')

for (i in seq(1,10)){

n <- 300 * i

set.seed(324)
A <- sample(seq(100,10000000), size = n, replace = TRUE)
set.seed(987)
B <- sample(seq(1000,100000000), size = n, replace = TRUE)

C = cbind(A, B)
D = kernel_fn(as.matrix(dist(C)))

t_1 <- proc.time()

inv_sigma_obs <- GPCERF:::compute_inverse(D)

t_2 <- proc.time()

df <- rbind(df, data.frame(size=n, wc=t_2[[3]] - t_1[[3]]))

}
#print(inv_sigma_obs)

print(df)

library(ggplot2)

ggplot(data = df) + geom_line(aes(size, wc), color="red") +
                    geom_point(aes(size, wc), color="blue")
