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

// [[Rcpp::depends(RcppEigen)]]
#ifndef TFIDF_
#define TFIDF_ 1

#include "matrix.hpp"

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

Vec tf_bool(    const Vec &tfs, const double &L );
Vec tf_raw(     const Vec &tfs, const double &L );
Vec tf_norm(    const Vec &tfs, const double &L );
Vec tf_logNorm( const Vec &tfs, const double &L );
Vec tf_05norm(  const Vec &tfs, const double &L );
Vec tf_Knorm(   const Vec &tfs, const double &L, const double &K );

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

Vec idf_unary(  const Vec &dfs, const double &D );
Vec idf_plain(  const Vec &dfs, const double &D );
Vec idf_smooth( const Vec &dfs, const double &D );
Vec idf_max(    const Vec &dfs, const double &D );
Vec idf_prob(   const Vec &dfs, const double &D );

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

#endif // TFIDF_
