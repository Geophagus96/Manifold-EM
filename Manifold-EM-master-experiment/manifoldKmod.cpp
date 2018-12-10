#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
/*The new modified EM algorithm is adjusted for the number of clusters*/

//[[Rcpp::plugins(cpp11)]]
/*Since my Rcpp package is an old version, it only supports cpp11. If you are using the latest Rcpp version please change the cpp plugin to cpp14*/
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
/*assigning all points to the initial point that is the nearest*/
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
/*mean estimation for the M-Step*/
int mean_est(const mat&dist, const uvec&num_seq){
  int k = dist.n_cols;
  vec ones(k,fill::ones);
  vec sumdist = dist*ones;
  uvec center_ind = find(sumdist == min(sumdist));
  int center = num_seq(center_ind(0));
  return center;
}

//[[Rcpp::export]]
/*Variance estimation for the M-Step*/
double variance_est(const mat&dist, const int&center,const uvec&num_seq){
  mat center_id = dist.row(center);
  mat center_dist = center_id.cols(num_seq);
  mat c_var = pow(center_dist,2);
  double sigma = mean(mean(c_var));
  return sigma;
}

//[[Rcpp::export]]
umat cats_KEM(const mat&dist, const uvec&initials, const int & manifolds, const int &categories, const int &max_iter){
  int m = categories;
  int n = dist.n_cols;
  umat result_label(n,max_iter);
  uvec cats = nnc(dist.rows(initials));
  vec proportions(categories, fill::ones);
  double cates = categories;
  proportions = (1/cates)*proportions;
  uvec centers = initials;
  vec sigmas(categories,fill::ones);
  uvec allcate = linspace<uvec>(0,(m-1),m);
  for (int j = 0; j<max_iter; j++){
    uvec centers(m);
    vec sigmas(m);
    cout << allcate << endl;
    for (int k =0; k<m; k++){
      uvec subsign = find(cats == allcate(k));
      mat subdist1 = dist.cols(subsign);
      mat subdist2 = subdist1.rows(subsign);
      int center = mean_est(subdist2,subsign);
      centers(k) = center;
      sigmas(k) = variance_est(dist,center,subsign);
    }
    allcate = linspace<uvec>(0,(m-1),m);
    mat subdist = dist.cols(centers);
    vec denominator = 1/pow(sigmas,(manifolds/2));
    mat likelihood = exp(-0.5*subdist*diagmat(denominator))*diagmat(denominator)*diagmat(10*proportions);
    for (int i=0; i<n; i++){
      uvec cat = find(likelihood.row(i) == max(likelihood.row(i)));
      cats(i) = cat(0);
    }
    vec temp_prop(m, fill::zeros);
    for (int h = 0 ;h < m; h++){
      uvec sumi = find(cats == h);
      temp_prop(h) = double(sumi.n_elem)/double(n);    
    }
    
    uvec valid_prop = find(temp_prop >= 1e-8);
    proportions = temp_prop(valid_prop);
    m = proportions.n_elem;
    uvec temp_cate = allcate(valid_prop);
    allcate = temp_cate;
    result_label.col(j) = cats;
  }  
  return result_label;
}
