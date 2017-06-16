// [[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include <stdio.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include "mattools.hpp"

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

typedef Eigen::MatrixXd Mat;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;
typedef Eigen::VectorXd Vec;

/*
 ' Compute a (P)PMI matrix from raw frequency counts, stored as a column-oriented sparse matrix.
 '
 ' @param m    A colun-oriented sparse matrix containing cooccurrence counts
 ' @param ppmi A logical value indicating whether negative values should be truncated to 0 (i.e.
 '             compute PPMI instead of PMI). TRUE by default.
 ' @param ow   A logical value indicating wether the result should be destructively copied over the
 '             input matrix. TRUE by default.
 '
 ' @return The PMI value for each cell is equal to log( p(i,j) / p(i)p(j) ), i.e. the log of the
 '         observed probability over the expected probability. The PPMI truncates negative values
 '         to 0. WARNING: In order to maintain sparsity, -Inf values in the non-truncated case are
 '         replaced by 0.
 '
*/
// [[Rcpp::export]]
S4 cooc_to_pmi( SEXP m_, bool ppmi = true, bool ow = true ) {
  SpMat m = as<MSpMat>( m_ );
  Vec rs = Vec::Zero( m.rows() );
  Vec cs = Vec::Zero( m.cols() );
  double N = m.sum();
  for( int i = 0; i < m.outerSize(); i++ ) {
    for( SpMat::InnerIterator it( m, i ); it; ++it ) {
      double v = it.value();
      rs[it.row()] += v / N;
      cs[it.col()] += v / N;
    }
  }
  SpMat tgt = ow ? m : as<MSpMat>( clone( m_ ) );
  for( int i = 0; i < m.outerSize(); i++ ) {
    for( SpMat::InnerIterator it( m, i ); it; ++it ) {
      int r = it.row();
      int c = it.col();
      double v = ( m.coeff( r, c ) / N ) / ( rs[r] * cs[c] ) + ( ppmi ? 1 : 0 );
      tgt.coeffRef( r, c ) = std::log<double>( v ).real();
    }
  }

  return Rcpp::wrap( tgt );
}

// [[Rcpp::export]]
NumericMatrix tfidf_weight( NumericMatrix m_, NumericVector df_, int mode = 1, bool ow = false ) {
  Mat m = as<Mat>( m_ );
  Vec df = as<Vec>( clone( df_ ) );
  if( m.rows() != df.rows() ) {

    // stop( std::sprintf( "Wrong dim for DF vector: got %d, wanted %d", df.rows(), m.cols() ) ;)
  }

  Vec L = m.colwise().sum();
  double D = df.sum();
  df = ( df / D ).array().log1p().matrix();

  Mat out = ow ? m : as<Mat>( clone( m_ ) );
  for( int i = 0; i < m.cols(); i++ ) {
    out.col( i ) /= L[i] * ( std::log1p( D / df[i] ) );
  }

  df = ( df / D ).array().log1p().matrix();
  Rcout << "out: (" << out.rows() << "," << out.cols() << ")" << std::endl;
  return wrap( out );
}

class IDF {
  public:
  double weight( long D, long d ) {
  }
};

/*** R
  m <- matrix( c( 1,1,1,1,2,0,0,2,1,0,0,3 ), nrow = 6, byrow = TRUE )
  rownames( m ) <- c( "this", "is", "a", "another", "sample", "example" )
  m
  df <- rowSums( m )
  tfidf_weight( m, df )
*/

