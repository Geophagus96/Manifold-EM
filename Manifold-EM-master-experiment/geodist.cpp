#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <array>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::plugins(cpp11)]]
/*Since my Rcpp package is an old version, it only supports cpp11. If you are using the latest Rcpp version please change the cpp plugin to cpp14*/
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
mat geodist(const mat & pathdist, const mat &l2dist, const mat & ccdist, const mat & ccrep, const uvec & ccs){
  mat geodist = pathdist;
  int n = pathdist.n_rows;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i < j){
        if (!isfinite(pathdist(i,j))){
          int cc1 = ccs(i);
          int cc2 = ccs(j);
          int rep1 = ccrep((cc1-1),(cc2-1));
          int rep2 = ccrep((cc2-1),(cc1-1));
          geodist(i,j) = l2dist((rep1-1),(rep2-1))+ pathdist(i,(rep1-1))+pathdist(j,(rep2-1));
        }
      }
      else if (i>j){
        geodist(i,j) = geodist(j,i);
      }
    }
  }
  return geodist;
}
/*
 The function for calculating geodesic distances when there are more than one connected components. The first 
 parameter is the shortest path matrix, the second parameter is the euclidean distance matrix and the third 
 parameter is the distance matrix between different connected components.
 
 The distance between two connected components i and j is defined as such:
                                      d(i,j) = min d(ci, cj), ci in i, cj in j
 
 Therefore the fourth parameter is the connected components representatives matrix, whose row number and column
 number is the same as the number of connected components. The (i,j) th entry of which is the node of the i th 
 component that is closest to the j th connected component, namely the ci chosen for the definition of the 
 distance d(i,j).
 
 The last parameter is a vector indicating which connected component each node is in.
 */