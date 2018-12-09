#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

/*
Selecting initial guess for EM algorithm by merging the important nodes
until there are exactly k left (we assume that there are k clusters)
*/
uvec IPwithmerge(const mat&knng, const mat&distm, const uvec&nnnums, const int&cats, const int&tresh){
  /*
  cats: number of clusters
  tresh: treshold for labeling a point as important, i.e. with many neighbors
  (number of neighbors)
  */

  int n = knng.n_rows;
  int k+ = knng.n_cols;
  uvec centers(cats, fill::zeros);

  uvec imnode_ind = find(nnnums >= tresh);
  uvec imnode_val = nnnums(imnode_ind);
  int imnode_num = imnodeval.n_elem;  
  mat tempdistm = {distm + (n+1) * eye(n,n), linspace<uvec>(0, (n - 1), n)};
  tempdistm = tempdistm(imnode_ind, {imnode_ind, n});
  
  for (int i = imnode_num; i > cats ; i--){
    uvec nodepair = closenp(tempdistm, i);
    tempdistm = mergeRD(tempdistm, nodepair, i);
    }

  centers = tempdistm.col(cats);
  return centers
}

uvec closenp(const mat&tempdistm, const int&s){
  /*
  To find nodes pair with closest distance
  every return a pair of nodes with smallest distance in the matrix
  input: tempdistm is a s-by-(s+1) matrix
  */
  uvec nodepair = {0, 0};
  mat tempd = tempdistm(linspace<uvec>(0, (s - 1), s), linspace<uvec>(0, (s - 1), s))
  int mindist = min(min(tempd));
  nodepair = find(tempd == mindist, 1, "first");
  return nodepair
}

mat mergeRD(const mat&tempdistm, const uvec&nodepair, const int&s){
  /*
  merge nodes randomly
  returns the new distance matrix with an extra column of indices
  */
  int del = 0;
  double r = randn(1,1);
  if (r >= 0.5){
    del = nodepair(0);
  }
  else{
    del = nodepair(1);
  }

  tempdistm.shed_row(del);
  tempdistm.shed_col(del);

  return tempdistm
}

mat mergeweighted(){
  /*
  merge nodes according to their weights
  returns the new distance matrix with an extra column of indices
  */
}