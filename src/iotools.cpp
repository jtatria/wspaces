#include "iotools.hpp"

using namespace Rcpp;

//' Load sparse matrices from disk.
//'
//' Reads a sparse matrix from binary file containing a sequence of (long, long, double) triplets.
//'
//' This function uses a memory mapped file to read the given binary file as an array of triplets
//' containing the coordinates and values for non-zero entries in a sparse matrix. This format is
//' compatible with the binary format used natively by Glove. The component values in the triplets
//' are used to construct the i, j and x vectors needed by Matrix::sparseMatrix.
//'
//' @param file A character vector indicating the path to the file to be read.
//'
//' @return A Matrix::sparseMatrix built from the triplets read from the given binary file.
//'
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
S4 load_spm( std::string file ) {
  SparseMatrix* spm = spm_read( file.c_str() );
  double* x = spm->entries;
  int*    i = spm->row_indices;
  int*    j = spm->col_indices;
  S4 out("dgTMatrix");
  IntegerVector iv( i, i + spm->entry_count );
  IntegerVector jv( j, j + spm->entry_count );
  NumericVector xv( x, x + spm->entry_count );
  return dgTMatrix( iv, jv, xv );
}

S4 dgTMatrix( IntegerVector i, IntegerVector j, NumericVector x ) {
  Function f = Environment::namespace_env("Matrix")["sparseMatrix"];
  return f( _["i"]=i, _["j"]=j, _["x"]=x, _["index1"]=false );
}

S4 dgCMatrix( IntegerVector i, IntegerVector p, NumericVector x ) {
  Function f = Environment::namespace_env("Matrix")["sparseMatrix"];
  return f( _["i"]=i, _["p"]=p, _["x"]=x, _["index1"]=false );
}
