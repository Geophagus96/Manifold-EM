#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

uvec sigpointsampling(const mat&knng, const int &t){
  int n = knng.n_rows;
  int k = knng.n_cols;
  uvec nnnums(n,fill::zeros);
  for (int i = 0; i<n; i++){
    for (int j = 0; j<k; j++){
      int m = knng(i,j);
      nnnums(m-1) =nnnums(m-1)+1;
    }
  }
  return nnnums;
  
}