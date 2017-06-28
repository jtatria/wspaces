// [[Rcpp::depends(Matrix)]]
#ifndef IOTOOLS_
#define IOTOOLS_ 1

#include <Rcpp.h>
extern "C" {
  #include "SparseMatrix.h"
}

Rcpp::S4 dgTMatrix( Rcpp::IntegerVector i, Rcpp::IntegerVector j, Rcpp::NumericVector x );
Rcpp::S4 dgCMatrix( Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::NumericVector x );

#endif // SPARSE_MATRIX_H_
