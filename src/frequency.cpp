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

/***************************************************************************************************
                            Term frequency weighting functions.
***************************************************************************************************/

#include "tools.h"

using namespace Rcpp;

const int TF_BOOLEAN = 0;
const int TF_RAW     = 1;
const int TF_NORM    = 2;
const int TF_LOGNORM = 3;
const int TF_05NORM  = 4;
const int TF_KNORM   = 5;

const int IDF_UNARY  = 0;
const int IDF_PLAIN  = 1;
const int IDF_SMOOTH = 2;
const int IDF_MAX    = 3;
const int IDF_PROB   = 4;

Vec tf_bool(    const Vec&, const double& );
Vec tf_raw(     const Vec&, const double& );
Vec tf_norm(    const Vec&, const double& );
Vec tf_logNorm( const Vec&, const double& );
Vec tf_05norm(  const Vec&, const double& );
Vec tf_Knorm(   const Vec&, const double&, const double& );

// TF dispatch function
inline Vec tf( int &mode, const Vec &tfs, const double &L ) {
    Vec (*func)( const Vec&, const double& );
    switch( mode ) {
        case TF_BOOLEAN : func = tf_bool;    break;
        case TF_RAW     : func = tf_raw;     break;
        case TF_NORM    : func = tf_norm;    break;
        case TF_LOGNORM : func = tf_logNorm; break;
        case TF_05NORM  : func = tf_05norm;  break;
        default: Rcpp::stop( "Unknown mode for TF function" );
    }
    return (*func)( tfs, L );
}

Vec idf_unary(  const Vec&, const double& );
Vec idf_plain(  const Vec&, const double& );
Vec idf_smooth( const Vec&, const double& );
Vec idf_max(    const Vec&, const double& );
Vec idf_prob(   const Vec&, const double& );

// IDF dispatch function
inline Vec idf( const int &mode, const Vec &dfs, const double &D ) {
    Vec (*func)( const Vec&, const double &D );
    switch( mode ) {
        case IDF_UNARY  : func = idf_unary; break;
        case IDF_PLAIN  : func = idf_plain; break;
        case IDF_SMOOTH : func = idf_smooth; break;
        case IDF_MAX    : func = idf_max; break;
        case IDF_PROB   : func = idf_prob; break;
        default: Rcpp::stop( "Unknown mode for IDF function" );
    }
    return (*func)( dfs, D );
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
RVecD tf( RVecD tfs, double L, int mode = 2 ) {
    return wrap( tf( mode, as<Vec>( tfs ), L ) );
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
RVecD idf( RVecD dfs, double D, int mode = 2 ) {
  return wrap( idf( mode, as<Vec>( dfs ), D ) );
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
RMatD tfidf( RMatD tf_, RVecD df_, int tf_mode = 2, int idf_mode = 2, bool ow = false ) {
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


// TF functions
inline Vec tf_bool( const Vec &tfs, const double &L ) {
    return tfs.cwiseSign();
}

inline Vec tf_raw( const Vec &tfs, const double &L ) {
    return tfs;
}

inline Vec tf_norm( const Vec &tfs, const double &L ) {
    return tfs / L;
}

inline Vec tf_logNorm( const Vec &tfs, const double &L ) {
    return ( ( tfs / L ).array() + 1 ).log().matrix();
}

inline Vec tf_05norm( const Vec &tfs, const double &L ) {
    return tf_Knorm( tfs, L, 0.5 );
}

inline Vec tf_Knorm( const Vec &tfs, const double &L, const double &K = 0.5 ) {
    return ( tfs / tfs.maxCoeff() ) * ( K + ( 1 - K ) );
}

// IDF functions
inline Vec idf_unary( const Vec &dfs, const double &D ) {
    return dfs.cwiseSign();
}

inline Vec idf_plain( const Vec &dfs, const double &D ) {
    return ( dfs.array().inverse() * D ).log().matrix();
}

inline Vec idf_smooth( const Vec &dfs, const double &D ) {
    return ( ( dfs.array().inverse() * D ) + 1 ).log().matrix();
}

inline Vec idf_max( const Vec &dfs, const double &D ) {
    return ( dfs.array().inverse() * dfs.maxCoeff() ).log().matrix();
}

inline Vec idf_prob( const Vec &dfs, const double &D ) {
    return ( ( dfs.array().inverse() ) * ( ( dfs.array() ) * -1 + D ) ).log().matrix();
}
