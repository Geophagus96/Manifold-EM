#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

/*
Selecting initial guess for EM algorithm by merging the important nodes
until there are exactly k left (we assume that there are k clusters)
*/
//[[Rcpp::export]]
uvec IPmaxdistmerge(const mat&knng, const mat&distm, const uvec&nnnums, const int&cats, const int&tresh){
  /*
  cats: number of clusters
  tresh: treshold for labeling a point as important, i.e. with more than a certain number of neighbors
  output a uvec, denoting the indices of the centers
  */
  int n = knng.n_rows;
  int del;
  double maxdist;    
  uvec nodeind;
  uvec nodeind2;
  vec distsum = sum(distm);

  uvec imnode_ind = find(nnnums >= tresh);
  int imnode_num = imnode_ind.n_elem;
  mat tempdistm(n, n);
  tempdistm = distm + (n+1) * eye(n,n);
  tempdistm = tempdistm(imnode_ind, imnode_ind);

  for (int i = imnode_num; i > cats ; i--){    
    maxdist = max(max(tempdistm));
    nodeind = find(tempdistm == maxdist, 1, "first");

    int nodelabel = nodeind(0) + 1;
    int nodencol = nodelabel / i;
    int nodenrow = nodelabel - i * nodencol;
    if (nodenrow == 0){
      nodencol = nodencol - 1;
      nodenrow = i;
    }
    nodenrow = nodenrow - 1; 
    /*
    coordinate in tempdistm, instead of dism!!
    */

    if (distsum(nodenrow) < distsum(nodencol)){
      nodeind2 = find(tempdistm.col(nodenrow) == min(tempdistm.col(nodenrow)));
      del = nodeind2(0);
    }
    else{
      nodeind2 = find(tempdistm.col(nodencol) == min(tempdistm.col(nodencol)));
      del = nodeind2(0);
    }

    tempdistm.shed_row(del);
    tempdistm.shed_col(del);
    imnode_ind = imnode_ind(find(imnode_ind != imnode_ind(del)));
    }

  return imnode_ind;
}
