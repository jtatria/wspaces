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
                    Implementation of weigth/extent on sparse matrices
***************************************************************************************************/

// [[Rcpp::depends( RcppParallel )]]
#include "tools.hpp"
#include <functional>
#include <cmath>
#include <RcppParallel.h>

enum Weight {
    Type  = 0,
    Token = 1,
    PMI   = 2,
    WPMI  = 3,
    TTest = 4,
    ZTest = 5,
    AllR  = 6,
};

namespace impl {
    // Not used for now.
    struct CoocWorker : public RcppParallel::Worker {
        SpMat src;
        SpMat tgt;
        F<double,ind,ind> func;

        CoocWorker( const SpMat &src, SpMat &tgt, const F<double,ind,ind>& func )
            : src( src ), tgt( tgt ), func( func ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( ind i = begin; i < end; ++i ) {
                for( SpInIt srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
                    double v = srcIt.value();
                    ind r    = srcIt.row();
                    ind c    = srcIt.col();
                    tgtIt.valueRef() = func( v, r, c );
                }
            }
        }
    };

    void weight_cooc( SpMat &src, SpMat &tgt, const F<double,ind,ind> func ) {
        // Eigen::initParallel();
        // CoocWorker wrkr( src, tgt, func );
        // RcppParallel::parallelFor( 0, src.outerSize(), wrkr );
        for( ind i = 0; i < src.rows(); ++i ) {
            for( SpInIt srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
                double v = srcIt.value();
                ind r    = srcIt.row();
                ind c    = srcIt.col();
                tgtIt.valueRef() = func( v, r, c );
            }
        }
    }
};

F<double,ind,ind> type(
    const double N, const Vec& rmrg, const Vec& cmrg, bool pos
) {
    return [N,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c | W ) > 0 ? 1 : 0
        // TODO check math?
        return ( v / N ) / cmrg( c ) > 0 ? 1 : 0;
    };
}

F<double,ind,ind> token( const double N, const Vec& rmrg, const Vec& cmrg, bool pos ) {
    return [N,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c | W )
        // TODO check math?
        return ( v / N ) / cmrg( c );
    };
}

// All logs are computed as log( x + 1 ) to avoid non-zero values for zero entries.
F<double,ind,ind> pmi( const double N, const Vec& rmrg, const Vec& cmrg, const bool pos ) {
    const int add = pos ? 1 : 0;
    return [N,&rmrg,&cmrg,add]( const double& v, const ind& r, const ind& c ) -> double {
        // log( P( c, W ) / P( c ) * P( W ) )
        double obs = ( v / N );
        double exp = ( rmrg( r ) * cmrg( c ) );
        return std::log<double>( ( obs / exp ) + add ).real();
    };
}

// All logs are computed as log( x + 1 ) to avoid non-zero values for zero entries.
F<double,ind,ind> wpmi( const double N, const Vec& rmrg, const Vec& cmrg, const bool pos ) {
    const int add = pos ? 1 : 0;
    return [N,&rmrg,&cmrg,add]( const double& v, const ind& r, const ind& c ) -> double {
        // log( P( c, W ) / P( c ) * P( W ) ) * P( c, W )
        double obs = ( v / N );
        double exp = ( rmrg( r ) * cmrg( c ) );
        return std::log<double>( ( obs / exp ) + add ).real() * obs;
    };
}

F<double,ind,ind> ttest( const double N, const Vec& rmrg, const Vec& cmrg, const bool pos ) {
    return [N,&rmrg,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c, W ) - P( c )*P( W ) / sqrt( P( c, W ) / N )
        double num = ( v / N ) - ( rmrg( r ) * cmrg( c ) );
        double den = std::sqrt<double>( ( v / N ) / N ).real();
        return num / den;
    };
}

// z-test weights are not 0 for for zero-entries in the input matrix and these are never computed.
F<double,ind,ind> ztest( const double N, const Vec& rmrg, const Vec& cmrg, const bool pos ) {
    return [N,&rmrg,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c, W ) - P( c )*P( W ) / sqrt( ( P( c ) * P( W ) / N ) )
        double num = ( v / N ) - ( rmrg( r ) * cmrg( c ) );
        double den = std::sqrt<double>( ( rmrg( r ) * cmrg( c ) ) / N ).real();
        return num / den;
    };
}

double L( double k, double n, double x ) {
    return std::pow<double>( x, k ) * std::pow<double>( 1 - x, n - k );
};

F<double,ind,ind> allr( const double N, const Vec& rmrg, const Vec& cmrg, const bool pos ) {
    return [N,&rmrg,&cmrg]( const double& v, const ind& r, const ind &c ) -> double {
        double num = L( v, rmrg( r ) * N, cmrg( c ) );
        double den = L( v, rmrg( r ) * N, v / ( rmrg( r ) * N ) );
        return -2 * std::log<double>( num / den ).real();
    };
}

inline F<double,ind,ind> get_func(
        Weight d,
        const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
    switch( d ) {
        case Type  : return  type( N, rmrg, cmrg, pos ); break;
        case Token : return token( N, rmrg, cmrg, pos ); break;
        case PMI   : return   pmi( N, rmrg, cmrg, pos ); break;
        case WPMI  : return  wpmi( N, rmrg, cmrg, pos ); break;
        case TTest : return ttest( N, rmrg, cmrg, pos ); break;
        case ZTest : return ztest( N, rmrg, cmrg, pos ); break;
        case AllR  : return  allr( N, rmrg, cmrg, pos ); break;
        default: Rcpp::stop( "Unknown function requested" );
    }
}

//' Weight the given cooccurrence matrix.
//'
//' Computes one of several weighting and extent functions on the non-zero entries of the given
//' (sparse) cooccurrence matrix and returns the transformed (sparse) matrix.
//'
//' \code{weight_cooc} is the actual workhorse function which can be used directly with its default
//' values to weight global cooccurrence matrices. \code{weight_sample} is a thin wrapper for
//' weighting corpus samples that uses marginals from a global cooccurrence matrix.
//'
//' Available weighting functions:
//'
//' \itemize{
//'     \item{0: Type weight: \eqn{P( c | W ) > 0 ? 1 : 0} }
//'     \item{1: Token weight: \eqn{P( c | W )} }
//'     \item{2: PMI: \eqn{ log( P( c, W ) / P( c ) * P( W ) ) } }
//'     \item{3: Weighted PMI: \eqn{ P( c, W ) * log( P( c, W ) / P( c ) * P( W ) ) } }
//'     \item{4: t-Test: \eqn{ P( c, W ) - P( c )*P( W ) / \sqrt{ P( c, W ) / N } } }
//'     \item{5: z-Test: \eqn{ P( c, W ) - P( c )*P( W ) / \sqrt{ ( P( c ) * P( W ) / N ) } } }
//'     \item{6: Log-likelihood approximation: See below.}
//' }
//'
//' Mode 6 uses [citation needed] 'allr' log-likelihood approximation:
//' \eqn{ -2 \frac{
//'     L( F( w, c ), F( w ), F( c ) )
//' }{
//'     L( F( w, c ), F( w ), \frac{F( W, c )}{F( W )} )
//' } } with \deqn{ L(k,n,x) = x^{k} * ( 1 - x )^{n-k} }
//'
//' @param m        A sparse matrix with raw cooccurrence counts or some function thereof.
//' @param rowm     An optional vector of length == nrow( m ) with population focal (i.e. row) term
//'                 probabilities. NULL by default, computed from m.
//' @param colm     An optional vector of length == ncol( m ) with population context (i.e. column)
//'                 term probabilities. NULL by default, computed from m.
//' @param positive Logical. Truncate negative values to zero if reasonable (i.e. use log( p + 1 )
//'                 internally). TRUE by default.
//' @param ow       Logical. Operate destructively on m by overwriting it with the result. FALSE by
//'                 default.
//' @param mode     Weighting function to apply on m. See details.
//'
//' @return a sparse matrix with similar structure to m with the results of the weighting function.
//'
//' @export
// [[Rcpp::export]]
Rcpp::S4 weight_cooc(
    Rcpp::S4 m,
    Rcpp::Nullable<RVecD> rowm = R_NilValue,
    Rcpp::Nullable<RVecD> colm = R_NilValue,
    bool positive = true, bool ow = false, int mode = 2 // plain pmi
) {
    SpMat src = Rcpp::as<MSpMat>( m );
    Vec rmrg = rowm.isNull() ? marg_prb( src, Margin::Row ) : Rcpp::as<Vec>( rowm.get() );
    Vec cmrg = colm.isNull() ? marg_prb( src, Margin::Col ) : Rcpp::as<Vec>( colm.get() );
    if( rmrg.size() != src.cols() ) Rcpp::stop( "Wrong dimension for row marginal vector" );
    if( cmrg.size() != src.cols() ) Rcpp::stop( "Wrong dimension for column marginal vector" );
    // rmrg = ( rmrg.array() / rmrg.sum() ).matrix().eval();
    // cmrg = ( cmrg.array() / cmrg.sum() ).matrix().eval();
    SpMat tgt = ow ? src : Rcpp::as<MSpMat>( clone( m ) );
    F<double,ind,ind> func = get_func( static_cast<Weight>( mode ), src.sum(), rmrg, cmrg, positive );
    impl::weight_cooc( src, tgt, func );
    Rcpp::S4 out = Rcpp::wrap( tgt );
    out.attr( "Dimnames" ) = m.attr( "Dimnames" );
    return out;
}
