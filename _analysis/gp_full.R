library(xgboost)

# Step 1: Generate synthetic data
sim.data <- generate_synthetic_data(sample_size = 500, gps_spec = 3)

# Step 2: Compute true curve
tru.curve <- sapply(seq(0,20,0.1), function(w){
  tru_R(w, sim.data)
})

# Step 3: Estimate GPS function
e_gps <- xgboost(label=sim.data$treat, data=as.matrix(sim.data[,-(1:2)]), nrounds = 50)
e_gps_pred <- predict(e_gps,as.matrix(sim.data[,-(1:2)]))
e_gps_std <- sd(sim.data$treat-e_gps_pred)
GPS <- dnorm(sim.data$treat, mean = e_gps_pred, sd = e_gps_std, log = T)


# Step 4: Get exposure values
w = 0.2 #seq(0,20,0.1)

# Step 5: Hyper parameters
alpha = 0.02
beta = 0.04
gamma_over_sigma = 1

# Step 6: Convert input data into matrix
x.design = model.matrix(~cf1+cf2+cf3+cf4+cf5+cf6-1, data = sim.data)

# Step 7: Define kernel function
kernel.fn = function(x) exp(-x^2)

# Step 8: Compute The inverse value (K + sigma^2 * I)^-1
treat <- sim.data$treat
obs.use = cbind(treat*sqrt(1/alpha), GPS*sqrt(1/beta))
Sigma.obs = gamma_over_sigma*kernel.fn(as.matrix(dist(obs.use))) + diag(nrow(obs.use)) # n*n
inv.Sigma.obs = chol2inv(chol(Sigma.obs)) # n*n

# Step 9: Compute GPS for new w
GPS.new = stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T) # Vector of length n

# Step 10: Compute Kappa
# Sigma.cross = kappa/sigma^2 : Is always n^2 matrix.
# each column of Sigma.cross is ki.
w.obs = sim.data$treat
obs.new = cbind(w*sqrt(1/alpha), GPS.new*sqrt(1/beta))
Sigma.cross = gamma_over_sigma*kernel.fn(spatstat.geom::crossdist(obs.new[,1],
                                                                  obs.new[,2],
                                                                  obs.use[,1],
                                                                  obs.use[,2])) # n*n

# Step 11: Compute weight
# weight is the same as a matrix in the paper
# weigts.final = invers of paranthesis * kappa
weights.final <- c((rep(1/length(w.obs),length(w.obs))%*%Sigma.cross)%*%inv.Sigma.obs)
weights.final[weights.final<0] = 0
weights.final = weights.final/sum(weights.final)



# Step 12: Compute est  (est is the same as m in the paper)
Y <- sim.data$Y
est = Y%*%weights.final

# Step 13: Compute rho
# this computes rho_r(w) for each covariate r
w.mean = sum(sim.data$treat*weights.final)
w.sd = sqrt(sum((sim.data$treat - w.mean)^2*weights.final))
w.stan = (sim.data$treat - w.mean)/w.sd

x.mean = colMeans(x.design*weights.final)
x.cov = (t(x.design) - x.mean)%*%diag(weights.final)%*%t(t(x.design) - x.mean)
x.stan = t(t(solve(chol(x.cov)))%*%(t(x.design) - x.mean))
col.all <- c(abs(c(t(x.stan)%*%diag(weights.final)%*%w.stan)), est )
tune_res <- list(cb = mean(col.all[1:6]), est = col.all[7])

# Step 14: return covariate balance and est.
rho = mean(tune_res$cb)
est =  tune_res$est
