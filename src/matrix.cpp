/*
 * Copyright (C) 2017 José Tomás Atria <jtatria at gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////// Matrix manipulation tools /////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include "matrix.hpp"
#include "tfidf.hpp"

using namespace Rcpp;

const double EPSILON = 0.00000001;

//' (Positive) Pointwise mutual information.
//'
//' Compute a (P)PMI matrix from raw frequency counts, stored as a column-oriented sparse matrix.
//'
//' The PMI value for each cell is equal to log( p(i,j) / p(i)p(j) ), i.e. the log of the observed
//' probability over the expected probability.
//'
//' The PPMI variant adds 1 to this value in order to truncate all values to 0 and maintain
//' sparsity.
//'
//' If rs or cs are given, the computation of the expected probabilities will use the values
//' contained in these vectors. This allows computing correct PMI values for partial cooccurrence
//' matrices, e.g. those built over a subsample of the entire corpus. If no marginal vectors are
//' given, these are computed from the row sums and column sums of the cooccurrence matrix,
//' respectively; this is equivalent to computing an 'empirical' PMI value from the given set of
//' observations, instead of computing it using the global corpus statistics.
//'
//' The computation of observed probabilities is obviously always based on the total number of
//' observations in the given cooccurrence matrix, independently of the marginal vectors used for
//' the expected probabilities.
//'
//' WARNING: In order to prevent a memory explosion, zero-values in the input matrix \emph{are never
//' calculated}, which has the numerical side-effect of replacing all -Inf values for the plain PMI
//' version with zeros. I have not yet thought of a better solution to distinguish -Infs from
//' actual 0s in the computed result.
//'
//' @param m     A colun-oriented sparse matrix containing cooccurrence counts.
//' @param rs    A numeric vector for row marginals (i.e. total DF for each term).
//' @param cs    A numeric vector for column marginals (i.e. total TF for each corpus segment).
//' @param ppmi  A logical value indicating whether negative values should be truncated to 0 (i.e.
//'              compute PPMI instead of PMI). TRUE by default.
//' @param ow    A logical value indicating wether the result should be destructively copied over
//'              the input matrix. FALSE by default.
//'
//' @return An (column-stored sparse) matrix, isomorphic to m_ with the (P)PMI values for m_.
//'         If ow == TRUE, m_ is replaced with this value.
//' @export
//' @importFrom Matrix rowSums colSums
// [[Rcpp::export]]
S4 cooc_to_pmi( S4 m_, NumericVector rs_, NumericVector cs_, bool ppmi = true, bool ow = false ) {
    SpMat src = as<MSpMat>( m_ );
    Vec   rs  = as<Vec>( rs_ );
    Vec   cs  = as<Vec>( cs_ );
    if( rs.size() != src.rows() ) {
        Rcpp::stop( "Wrong dimensions for row marginal vector" );
    }
    if( cs.size() != src.cols() ) {
        Rcpp::stop( "Wrong dimensions for col marginal vector" );
    }
    if( rs.sum() - 1 > EPSILON || cs.sum() - 1 > EPSILON ) {
        Rcpp::stop( "Marginals don't look like probability disributions" );
    }
    double N  = src.sum();
    SpMat tgt = ow ? src : as<MSpMat>( clone( m_ ) );
    for( int i = 0; i < src.outerSize(); i++ ) {
        for( SpMat::InnerIterator srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
            int r = srcIt.row();
            int c = srcIt.col();
            // Compiler abuse to facilitate debuging. These should all be optimized out.
            double obs = ( srcIt.value() / N );
            double exp = rs[r] * cs[c];
            double rat = obs / exp;
            double v   = std::log<double>( rat + ( ppmi ? 1 : 0 ) ).real();
            tgtIt.valueRef() = v;
        }
    }
    S4 out = wrap( tgt );
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
    Mat m  = as<Mat>( tf_ );
    Vec df = as<Vec>( clone( df_ ) );
    Vec L  = m.colwise().sum();
    Vec idfs = idf( idf_mode, df, m.cols() );
    Mat out = ow ? m : as<Mat>( clone( tf_ ) );
    for( int i = 0; i < m.cols(); i++ ) {
        out.col( i ) = ( tf( tf_mode, m.col( i ), L[i] ) ).array() * ( idfs.array() );
    }
    return wrap( out );
}

//' Compute IDF weights for the given DF vector.
//'
//' Compute IDF weigghts for the given DF vector over the given D document sample size, using the
//' given weighting mode.
//'
//' \itemize{
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
//' @param dfs  A vector with raw document frequencies.
//' @param D    Numeric vector of length one indicating the document sample size.
//' @param mode Numeric vector of length one indicating the mode used for IDF calculation. See
//'             details.
//' @return A vector of the same length as dfs with IDF weights.
//'
// [[Rcpp::export]]
NumericVector idf( NumericVector dfs, double D, int mode = 2 ) {
  return wrap( idf( mode, as<Vec>( dfs ), D ) );
}

//' Compute TF weights for the given raw TF vector.
//'
//' Compute TF weights for the given raw TF vector over the given L document length, using the
//' given weighting mode.
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
//' }
//'
//' @param tfs  A vector with raw term frequencies.
//' @param L    Numeric vector of length one, indicating the document's total term length.
//' @param mode Numeric vector of length one indicating the mode used for TF calculation. See
//'             details.
//' @return A vector of the same size as tfs with TF weights.
//'
// [[Rcpp::export]]
NumericVector tf( NumericVector tfs, double L, int mode = 2 ) {
    return wrap( tf( mode, as<Vec>( tfs ), L ) );
}

//' Cosine distances.
//'
//' Compute the cosine distance between the given vectors. The cosine distance between vectors
//' i and j is equal to the dot product between them divided by the product of their norms:
//' \eqn{ ( i \dot j ) / ||i||*||j|| }.
//'
//' The given matrix will be interpreted as an array of column vectors for which cosine distances
//' will be computed. If transpose == TRUE, the matrix will be interpreted as an array of row
//' vectors. In any case, the resulting matrix will be square, symmetric, and with dimensions
//' equal to the number of given (column or row) vectors.
//'
//' @param m_        A matrix of column (or row) vectors for which to compute distances.
//' @param transpose A logical vector indicating if distances should be computed across row vectors.
//'                  False by default.
//'
//' @return A square symmetric matrix of dimension equal to the column (or row) dimension of the
//'         given input matrix of cosine distances between the column (or row) vectors in the input
//'         matrix.
// [[Rcpp::export]]
NumericMatrix vector_cosine( NumericMatrix m_, bool transpose = false ) {
    Mat m = as<Mat>( m_ );
    m = transpose ? m.transpose() : m;
    Mat out( m.cols(), m.cols() );
    for( int i = 0; i < m.cols(); i++ ) {
        for( int j = 0; j < m.cols(); j++ ) {
            Vec vi = m.col( i );
            Vec vj = m.col( j );
            out.coeffRef( i, j ) = ( vi.dot( vj ) ) / ( vi.norm() * vj.norm() );
        }
    }
    return wrap( out );
}

/*** R
    m <- matrix( c( 1,1,1,1,2,0,0,2,1,0,0,3 ), nrow = 6, byrow = TRUE )
    rownames( m ) <- c( "this", "is", "a", "another", "sample", "example" )
    m
    df <- rowSums( ( m > 0 ) * 1 )
    tfidf_weight( m, df )
*/
