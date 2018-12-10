#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::plugins(cpp11)]]
/*Since my Rcpp package is an old version, it only supports cpp11. If you are using the latest Rcpp version please change the cpp plugin to cpp14*/
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

/*
 The function of calculating the distance between different connected components.
 */
uvec minccdist(const mat & ccdist){
  uvec x = {0,0};
  double minmum = ccdist(0,0);
  int k1 = ccdist.n_rows;
  int k2 = ccdist.n_cols;
  for (int i = 0; i < k1; i++){
    for (int j = 0; j < k2; j++){
      if (ccdist(i,j) < minmum){
        minmum = ccdist(i,j);
        x[0] = i;
        x[1] = j;
      }
    }
  }
  return x;
}
