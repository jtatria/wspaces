// [[Rcpp::depends(Matrix)]]
#include <Rcpp.h>
#include "iotools.hpp"

extern "C" {
  #include "SparseMatrix.h"
}

using namespace Rcpp;

//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
S4 load_spm( CharacterVector file ) {
  SparseMatrix* spm = spm_read( file[0] );
  double* x = spm->entries;
  int*    i = spm->row_indices;
  int*    j = spm->col_indices;
  S4 out("dgTMatrix");
  IntegerVector iv( i, i + spm->entry_count );
  IntegerVector jv( j, j + spm->entry_count );
  NumericVector xv( x, x + spm->entry_count );
  return build_dgTMatrix( iv, jv, xv );
}

SEXP build_dgTMatrix( IntegerVector i, IntegerVector j, NumericVector x ) {
  Function f = Environment::namespace_env("Matrix")["sparseMatrix"];
  return f( _["i"]=i, _["j"]=j, _["x"]=x, _["index1"]=false );
}
