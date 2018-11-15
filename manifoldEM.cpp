#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
uvec nnc(const mat&dist){
  int n = dist.n_cols;
  int k = dist.n_rows;
  uvec cate(n);
  for (int i = 0; i <n; i++){
    vec disti = dist.col(i);
    int c = 0;
    for (int j = 1; j < k; j++){
      if (disti(j)<disti(c)){
        c = j;
        }
      cate(i) = c;
      }
  }
  return cate;
}

//[[Rcpp::export]]
uvec mean_est(const mat&dist, const uvec&num_seq){
  int k = dist.n_cols;
  vec ones(k,fill::ones);
  vec sumdist = dist*ones;
  uvec center_ind = find(sumdist == min(sumdist));
  uvec center = num_seq(center_ind);
  return center;
}

//[[Rcpp::export]]
double variance_est(const mat&dist, const uvec&center,const uvec&num_seq){
  mat center_id = dist.row(center(0));
  mat center_dist = center_id.cols(num_seq);
  mat c_var = pow(center_dist,2);
  double sigma = mean(mean(c_var));
  return sigma;
}


//[[Rcpp::export]]
uvec cats_EM(const mat&dist, const uvec&initials, const int & manifolds, const int &categories, const int &max_iter){
  int n = dist.n_cols;
  uvec cats = nnc(dist.rows(initials));
  vec proportions(categories, fill::ones);
  /*proportions = 0.25*proportions;*/
  uvec centers = initials;
  vec sigmas(categories,fill::ones);
  for (int j = 0; j<max_iter; j++){
    for (int k =0; k<categories; k++){
      uvec subsign = find(cats == k);
      mat subdist1 = dist.cols(subsign);
      mat subdist2 = subdist1.rows(subsign);
      uvec center = mean_est(subdist2,subsign);
      centers(k) = center(0);
      sigmas(k) = variance_est(dist,center,subsign);
    }
    mat subdist = dist.cols(centers);
    vec denominator = 1/pow(sigmas,(manifolds/2));
    mat likelihood = exp(-0.5*subdist*diagmat(denominator))*diagmat(denominator)/*diagmat(10*proportions)*/;
    for (int i=0; i<n; i++){
      uvec cat = find(likelihood.row(i) == max(likelihood.row(i)));
      cats(i) = cat(0);
    }
    /*for (int h = 0 ;h < categories; h++){
      uvec sumi = find(cats == h);
      proportions(h) = sumi.n_elem/double(n);    }
  }*/
  }  
  return cats;
}  