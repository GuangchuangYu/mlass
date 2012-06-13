#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP computeCentroids(SEXP X, SEXP idx, SEXP K) {
  NumericMatrix Xm(X);
  NumericVector ridx(idx);
  int rK = as<int>(K);

  int nrows = Xm.nrow();
  int ncols = Xm.ncol();

  NumericMatrix centroids(rK, ncols);
  NumericVector cnt(rK);

  for (int j=0; j < ncols; j++) {
    for (int k=0; k < rK; k++) {
      for (int i=0; i < nrows; i++) {
        if (ridx(i) == k+1) {
          centroids(k,j) += Xm(i,j);
          cnt(k) +=1;
        }
      }
    }

  }

  for (int k=0; k < rK; k++) {
    cnt(k) /= ncols;
    for (int j=0; j < ncols; j++) {
      centroids(k,j) /= cnt(k);
    }
  }
  return Rcpp::wrap(centroids);
}
