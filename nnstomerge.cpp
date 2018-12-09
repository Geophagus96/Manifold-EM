#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

int closenp(const mat&tempdistm);

/*
Selecting initial guess for EM algorithm by merging the important nodes
until there are exactly k left (we assume that there are k clusters)
*/
//[[Rcpp::export]]
uvec IPstomerge(const mat&knng, const mat&distm, const uvec&nnnums, const int&cats, const int&tresh){
  /*
  cats: number of clusters
  tresh: treshold for labeling a point as important, i.e. with more than a certain number of neighbors
  output a uvec, denoting the indices of the centers
  */
  int n = knng.n_rows;
  int del;
  double r;

  uvec imnode_ind = find(nnnums >= tresh);
  int imnode_num = imnode_ind.n_elem;
  mat tempdistm(n, n);
  tempdistm = distm + (n+1) * eye(n,n);
  tempdistm = tempdistm(imnode_ind, imnode_ind);

  for (int i = imnode_num; i > cats ; i--){
    int nodelabel = closenp(tempdistm) + 1;
    int nodencol = nodelabel / i;
    int nodenrow = nodelabel - i * nodencol;
    if (nodenrow == 0){
      nodencol = nodencol - 1;
      nodenrow = i;
    }    
    r = randu();
    if (r >= 0.5){
      del = nodenrow - 1;
    }
    else{
      del = nodencol;
    }
    tempdistm.shed_row(del);
    tempdistm.shed_col(del);
    imnode_ind = imnode_ind(find(imnode_ind != imnode_ind(del)));
    }

  return imnode_ind;
}

int closenp(const mat&tempdistm){
  /*
  To find nodes pair with closest distance
  every return a pair of nodes with smallest distance in the matrix
  input: tempdistm is a s-by-s matrix
  */
  uvec nodepair;
  double mindist;
  mindist = min(min(tempdistm));
  nodepair = find(tempdistm == mindist, 1, "first");
  int out = nodepair[0];
  return out;
}
