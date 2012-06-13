#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP findClosestCentroids(SEXP X, SEXP centroids) {
  NumericMatrix Xm(X);
  NumericMatrix center(centroids);
  
  int K = center.nrow();
  int nrows = Xm.nrow();
  int ncols = Xm.ncol();
  
  std::vector<int> idx(nrows);
  
  for (int i=0; i < nrows; i++) {
    
    int min_idx;
    float min;

    for (int k=0; k < K; k++) {
      float s = 0.0;  
      for (int j=0; j < ncols; j++) {
        s += pow( Xm(i,j) - center(k,j) , 2 );
     }

      if (k == 0) {
        min = s;
        min_idx = k+1;
      } else if (min > s) {
         min = s;
         min_idx = k+1;
       }

    }

    idx[i] = min_idx;
  }

  return Rcpp::wrap(idx);
}
