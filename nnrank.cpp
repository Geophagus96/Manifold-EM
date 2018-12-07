#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]
/*Rank all data points according to how many neighbourhoods they are in*/
//[[Rcpp::export]]
uvec sigpointsampling(const mat&knng){
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

/*Selecting all initial data points using the NN-rank method*/
//[[Rcpp::export]]
uvec initialpointsamp(const mat&knng, const uvec&nnnums, const int &cats){
  uvec nnnum_id= sort_index(nnnums,1);
  int n = knng.n_rows;
  int k = knng.n_cols;
  int center_num = 0;
  uvec centers(cats,fill::zeros);
  centers(0) = nnnum_id(0);
  center_num = center_num + 1;
  for (int i = 1; i< n; i++){
    int temp_id = nnnum_id(i);
    int new_c = 1;
    for (int j = 0;j< center_num; j++){
      for (int l = 0; l<k; l++){
        if (knng(centers(j),l) == temp_id){
          new_c = 0;
          break;
          break;
        }
      }
    }
    if (new_c == 1){
      centers(center_num) = temp_id;
      center_num = center_num+1;
    }
    if (center_num == cats){
      break;
    }
  }
  return centers;
}
