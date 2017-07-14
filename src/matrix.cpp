#include <math.h>
#include <stdio.h>
#include "tfidf.hpp"

using namespace Rcpp;

//' @internal
// [[Rcpp::export]]
S4 spm_pmi( S4 m_, NumericVector rs_, NumericVector cs_, bool ppmi = true, bool ow = false ) {
  SpMat src = Rcpp::as<MSpMat>( m_ );
  Vec   rs  = Rcpp::as<Vec>( rs_ );
  Vec   cs  = Rcpp::as<Vec>( cs_ );

  double N = src.sum();
  SpMat tgt = ow ? src : Rcpp::as<MSpMat>( Rcpp::clone( m_ ) );
  for( int i = 0; i < src.outerSize(); i++ ) {
    for( SpMat::InnerIterator srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
      int r = srcIt.row();
      int c = srcIt.col();
      tgtIt.valueRef() = std::log<double>(
        ( srcIt.value() / N ) / ( rs[r] * cs[c] ) + ( ppmi ? 1 : 0 )
      ).real();
    }
  }

  S4 out = Rcpp::wrap( tgt );
  out.slot( "Dimnames" ) = m_.slot( "Dimnames" );
  return out;
}

//' TF-IDF weighting.
//'
//' Weigth the given frequency matrix using the TF-IDF stategy.
//'
//' TF-IDF weights attempt to moderate the effect of very common terms, by dividing the total
//' frequency of a term within a context (i.e. a document, but in general any corpus segment),
//' by the number of contexts in which the term appears.
//'
//' The idea behind this approach is that terms that appear in every possible context do not
//' provide any additional information to the contexts in which they appear. Hence, all TF-IDF
//' strategies compute weights as some variation of \eqn{TF_{t,c} / IDF_{t}}, where
//' \eqn{TF_{t,c}} is a monotonic function of a term t's prevalence within a specific context
//' c and \eqn{IDF_{t}} is an inversely monotonic function of the term t's prevalence in all
//' contexts across the entire corpus.
//'
//' Note that the TF term is valid for a term in a context, while the IDF term is valid for a
//' term across the entire corpus.
//'
//' Values of tf_mode and idf_mode indicate how the TF and IDF terms are computed, as indicated
//' below.
//'
//' \itemize{
//'   \item{TF modes}
//'   \itemize{
//'     \item{0: Boolean: 1 if tf > 0; 0 otherwise.}
//'     \item{1: Raw: Raw term frequency.}
//'     \item{2: Normalized: (default) Term frequency divided by the total number of terms
//'              in document (the 'length').}
//'     \item{3: Log-normlized: Natural log of the term frequency over total document
//'              terms, + 1.}
//'     \item{4: 0.5 normalized: K*(1-K) * (tf / max( tf ) ), with K set to 0.5.}
//'   }
//'   \item{IDF modes}
//'   \itemize{
//'     \item{0: Unary: 1 if df > 0, but terms with 0 DF are by definition excluded of the
//'              lexicon, so 1.}
//'     \item{1: Plain: log of total number of documents, D, over the term's df.}
//'     \item{2: Smooth: (default) log of total number of documents, D, over the term's df,
//'              plus 1.}
//'     \item{3: Max: log of maximum df, over the term's df.}
//'     \item{4: Probabilistic: log of the total number of documents minus the term's df
//'              over the term's df.}
//'   }
//' }
//'
//' @param tf_       A matrix with one row for each term, and as many columns as documents or corpus
//'                  segments there are frequencies for.
//' @param df_       A vector of length equal the number of rows in tf_, containing document
//'                  frequencies, to compute the IDF component.
//' @param tf_mode   A term frequency weigthing strategy. See details.
//' @param idf_mode  An inverse document frequency weigthing strategy. See details.
//' @param ow        A logical vector indicating whether the result should be destructively copied
//'                  over the input matrix.
//'
//' @return An isomorphic matrix to tf_, with entries weighted by the given strategy. If ow == TRUE,
//'        tf_ is replaced with this value.
// [[Rcpp::export]]
NumericMatrix weight_tfidf(
    NumericMatrix tf_, NumericVector df_, int tf_mode = 2, int idf_mode = 2, bool ow = false
) {
  if( tf_.rows() != df_.length() ) {
    // stop( std::sprintf( "Wrong dim for DF vector: got %d, wanted %d", df.rows(), m.cols() ) ;)
  }
  Mat m  = Rcpp::as<Mat>( tf_ );
  Vec df = Rcpp::as<Vec>( clone( df_ ) );

  Vec L  = m.colwise().sum();
  Vec idfs = idf( idf_mode, df, m.cols() );
  Mat out = ow ? m : Rcpp::as<Mat>( Rcpp::clone( tf_ ) );
  for( int i = 0; i < m.cols(); i++ ) {
    out.col( i ) = ( tf( tf_mode, m.col( i ), L[i] ) ).array() * ( idfs.array() );
  }

  return Rcpp::wrap( out );
}

// [[Rcpp::export]]
Rcpp::NumericVector idf( Rcpp::NumericVector dfs_, double D, int mode = 2 ) {
  return Rcpp::wrap( idf( mode, Rcpp::as<Vec>( dfs_ ), D ) );
}

// [[Rcpp::export]]
Rcpp::NumericVector tf( Rcpp::NumericVector tfs_, double L, int mode = 2 ) {
  return Rcpp::wrap( tf( mode, Rcpp::as<Vec>( tfs_ ), L ) );
}

//' Cosine distances.
//'
//' Compute the cosine distance between the given vectors. The cosine distance between vectors
//' i and j is equal to the dot product between them divided by the product of their norms:
//' \eqn{ ( i \dot j ) / [[i]]*[[j]] }.
//'
//' The given matrix will be interpreted as an array of column vectors for which cosine distances
//' will be computed. If transpose == TRUE, the matrix will be interpreted as an array of row
//' vectors. In any case, the resulting matrix will be square, symmetric, and with dimensions
//' equal to the number of (column or row) vectors.
//'
//'
//' @param m_        A matrix of column (or row) vectors for which to compute distances.
//' @param transpose A logical vector indicating if distances should be computed across row vectors.
//'                  False by default.
//'
//' @return A square symmetric matrix of dimension equal to the column (or row) dimension of the
//'         given input matrix of cosine distances between the column (or row) vectors in the input
//'         matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix vector_cosine( Rcpp::NumericMatrix m_, bool transpose = false ) {
  Mat m = Rcpp::as<Mat>( m_ );
  m = transpose ? m.transpose() : m;
  Mat out( m.cols(), m.cols() );
  for( int i = 0; i < m.cols(); i++ ) {
    for( int j = 0; j < m.cols(); j++ ) {
      Vec vi = m.col( i );
      Vec vj = m.col( j );
      out.coeffRef( i, j ) = ( vi.dot( vj ) ) / ( vi.norm() * vj.norm() );
    }
  }
  return Rcpp::wrap( out );
}

/*** R
  m <- matrix( c( 1,1,1,1,2,0,0,2,1,0,0,3 ), nrow = 6, byrow = TRUE )
  rownames( m ) <- c( "this", "is", "a", "another", "sample", "example" )
  m
  df <- rowSums( ( m > 0 ) * 1 )
  tfidf_weight( m, df )
*/
