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

#include "tools.hpp"

using namespace Rcpp;

enum TF {
    Boolean = 0,
    Raw     = 1,
    Norm    = 2,
    Lognorm = 3,
    K5Norm  = 4,
};

enum IDF {
    Unary  = 0,
    Plain  = 1,
    Smooth = 2,
    Max    = 3,
    Prob   = 4,
};

// TF dispatch function
inline std::function<Vec(Vec,scalar)> get_tf( const int mode ) {
    TF tf = static_cast<TF>( mode );
    switch( tf ) {
        case Boolean : return []( const Vec& tf, const scalar L ) -> Vec {
            return tf.cwiseSign();
        };
        case Raw     : return []( const Vec& tf, const scalar L ) -> Vec {
            return tf;
        };
        case Norm    : return []( const Vec& tf, const scalar L ) -> Vec {
            return tf / L;
        };
        case Lognorm : return []( const Vec& tf, const scalar L ) -> Vec {
            return ( ( tf / L ).array() + 1 ).log().matrix();
        };
        case K5Norm  : return []( const Vec& tf, const scalar L ) -> Vec {
            scalar K = 0.5;
            return ( tf / tf.maxCoeff() ) * ( K + ( 1 - K ) );
        };
        default : Rcpp::stop( "Unknown mode for TF function" );
    }
}

// IDF dispatch function
inline std::function<Vec(Vec,scalar)> get_idf( const int mode ) {
    IDF idf = static_cast<IDF>( mode );
    switch( idf ) {
        case Unary  : return []( const Vec& df, const scalar D ) -> Vec {
            return df.cwiseSign();
        };
        case Plain  : return []( const Vec& df, const scalar D ) -> Vec {
            return ( df.array().inverse() * D ).log().matrix();
        };
        case Smooth : return []( const Vec& df, const scalar D ) -> Vec {
            return ( ( df.array().inverse() * D ) + 1 ).log().matrix();
        };
        case Max    : return []( const Vec& df, const scalar D ) -> Vec {
            return ( df.array().inverse() * df.maxCoeff() ).log().matrix();
        };
        case Prob   : return []( const Vec& df, const scalar D ) -> Vec {
            return ( ( df.array().inverse() ) * ( ( df.array() ) * -1 + D ) ).log().matrix();
        };
        default: Rcpp::stop( "Unknown mode for IDF function" );
    }
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
//' @export
// [[Rcpp::export]]
RVecD weight_tf( RVecD tfs, double L, int mode = 2 ) {
    return wrap( get_tf( static_cast<TF>( mode ) )( as<Vec>( tfs ), L ) );
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
//' @export
// [[Rcpp::export]]
RVecD weight_idf( RVecD dfs, double D, int mode = 2 ) {
  return wrap( get_idf( static_cast<IDF>( mode ) )( as<Vec>( dfs ), D ) );
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
//'              in the document or segment (the 'word length').}
//'     \item{3: Log-normlized: Natural log of the term frequency divided by the total number of
//'              terms in the document or segment (the 'word length'), + 1.}
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
//' @param tf       A matrix with one row for each term, and as many columns as documents or corpus
//'                 segments there are frequencies for.
//' @param fd       A vector of length equal the number of rows in tf_, containing document
//'                 frequencies, to compute the IDF component.
//' @param D        Total number of documents. Defaults to 0, in which case it is taken from the
//'                 number of columns in tf.
//' @param tf_mode  A term frequency weigthing strategy. See details.
//' @param idf_mode An inverse document frequency weigthing strategy. See details.
//' @param normal   Normalize the resulting vector to the (0,1] range.
//' @param ow       A logical vector indicating whether the result should be destructively copied
//'                 over the input matrix.
//'
//' @return An isomorphic matrix to tf_, with entries weighted by the given strategy. If ow == TRUE,
//'        tf_ is replaced with this value.
//' @export
// [[Rcpp::export]]
RMatD tfidf(
    RMatD tf, RVecD df, scalar D=0.0, int tf_mode = 2, int idf_mode = 2, bool normal=false,
    bool ow = false
) {
    if( tf.rows() != df.length() ) {
        // stop( std::sprintf( "Wrong dim for DF vector: got %d, wanted %d", df.rows(), m.cols() ) ;)
    }
    Mat src   = as<Mat>( tf );
    Mat tgt   = ow ? src : Mat( src.rows(), src.cols() );
    Vec dfs   = as<Vec>( df );
    Vec L     = src.colwise().sum();
    D = ( D <= 0.0 ? src.cols() : D );

    Vec idfs  = get_idf( idf_mode )( dfs, D );
    for( int i = 0; i < src.cols(); i++ ) {
        Vec res = get_tf( tf_mode )( src.col( i ), L[i] ).cwiseProduct( idfs );
        res = normal ? res / res.maxCoeff() : res;
        tgt.col( i ) = res;
    }

    RMatD out = wrap( tgt );
    out.attr( "Dimnames" ) = tf.attr( "Dimnames" );
    return wrap( out );
}

