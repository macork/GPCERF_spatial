#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat calc_cross(mat cross, mat within) {
  int n=cross.n_rows;
  mat sum(1,1,fill::zeros);
  for(int i=0;i<n;i++){
    sum += (cross.row(i)*within)*cross.row(i).t();
  }
  return(sum);
}