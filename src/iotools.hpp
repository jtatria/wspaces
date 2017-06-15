#ifndef IOTOOLS_
#define IOTOOLS_ 1

#include <Rcpp.h>

SEXP build_dgTMatrix( Rcpp::IntegerVector i, Rcpp::IntegerVector j, Rcpp::NumericVector x );

#endif // SPARSE_MATRIX_H_
