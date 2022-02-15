calc.ac = function(w, X, weights){
  w.mean = sum(w*weights)
  w.sd = sqrt(sum((w-w.mean)^2*weights))
  w.trans = (w-w.mean)/w.sd

  X.mean = colSums(X*weights)
  X.cov = (t(X) - X.mean)%*%diag(weights)%*%t(t(X)-X.mean)
  X.trans = t(t(solve(chol(X.cov)))%*%(t(X)-X.mean))

  c(w.trans%*%diag(weights)%*%X.trans)
}
