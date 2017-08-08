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
                Implementation of weigth/extent functions of cooccurrence counts.
***************************************************************************************************/


// [[Rcpp::depends( RcppParallel )]]
#include "tools.h"
#include <functional>
#include <cmath>
#include <RcppParallel.h>

using namespace Rcpp;

enum Weight {
    Type  = 0,
    Token = 1,
    PMI   = 2,
    WPMI  = 3,
    TTest = 4,
    ZTest = 5,
    AllR  = 6,
};

struct CoocWorker : public RcppParallel::Worker {
    SpMat src;
    SpMat tgt;
    F<double,ind,ind> func;

    CoocWorker( const SpMat& src, SpMat& tgt, const F<double,ind,ind>& func )
        : src( src ), tgt( tgt ), func( func ) {}

    void operator()( std::size_t begin, std::size_t end ) {
        for( ind i = begin; i < end; ++i ) {
            for( SpInIt srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
                double v = srcIt.value();
                ind r    = srcIt.row();
                ind c    = srcIt.col();
                // tgtIt.valueRef() = func( v, r, c );
                tgtIt.valueRef() = -1;
            }
        }
    }
};

F<double,ind,ind> type(
    const double N, const Vec& rmrg, const Vec& cmrg, bool pos
) {
    return [&,N,cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c | W ) > 0 ? 1 : 0
        return ( v / N ) / cmrg( c ) > 0 ? 1 : 0;
    };
}

F<double,ind,ind> token(
    const double N, const Vec& rmrg, const Vec& cmrg, bool pos
) {
    return [N,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c | W )
        return ( v / N ) / cmrg( c );
    };
}

// All logs are computed as log( x + 1 ) to avoid non-zero values for zero entries.
F<double,ind,ind> pmi(
    const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
    const int add = pos ? 1 : 0;
    return [N,&rmrg,&cmrg,add]( const double& v, const ind& r, const ind& c ) -> double {
        // log( P( c, W ) / P( c ) * P( W ) )
        double obs = ( v / N );
        double exp = ( rmrg( r ) * cmrg( c ) );
        return std::log<double>( ( obs / exp ) + add ).real();
    };
}

// All logs are computed as log( x + 1 ) to avoid non-zero values for zero entries.
F<double,ind,ind> wpmi(
    const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
    const int add = pos ? 1 : 0;
    return [N,&rmrg,&cmrg,add]( const double& v, const ind& r, const ind& c ) -> double {
        // log( P( c, W ) / P( c ) * P( W ) ) * P( c, W )
        double obs = ( v / N );
        double exp = ( rmrg( r ) * cmrg( c ) );
        return std::log<double>( ( obs / exp ) + add ).real() * obs;
    };
}

F<double,ind,ind> ttest(
    const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
    return [N,&rmrg,&cmrg]( const double& v, const ind& r, const ind& c ) -> double {
        // P( c, W ) - P( c )*P( W ) / sqrt( P( c, W ) / N )
        double num = ( v / N ) - ( rmrg( r ) * cmrg( c ) );
        double den = std::sqrt<double>( ( v / N ) / N ).real();
        return num / den;
    };
}

// z-test weights are not 0 for for zero-entries in the input matrix and these are never computed.
F<double,ind,ind> ztest(
    const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
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

F<double,ind,ind> allr(
    const double N, const Vec& rmrg, const Vec& cmrg, const bool pos
) {
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
        default: stop( "" );
    }
}

// [[Rcpp::export]]
S4 weight_cooc(
    S4 m, Nullable<RVecD> rowm = R_NilValue, Nullable<RVecD> colm = R_NilValue,
    bool positive = true, bool ow = false, int mode = 3 // plain pmi
) {
    SpMat src = as<MSpMat>( m );
    Vec rmrg = rowm.isNull() ? marg_prb( src, Margin::Row ) : as<Vec>( rowm.get() );
    Vec cmrg = colm.isNull() ? marg_prb( src, Margin::Col ) : as<Vec>( colm.get() );
    if( cmrg.size() != src.cols() ) stop( "Wrong dimension for row marginal vector" );
    if( cmrg.size() != src.cols() ) stop( "Wrong dimension for column marginal vector" );
    if( rmrg.sum() - 1 > EPSILON ) warning( "Row marginal doesn't seem to be a pr. distribution" );
    if( cmrg.sum() - 1 > EPSILON ) warning( "Col marginal doesn't seem to be a pr. distribution" );
    double N  = src.sum();
    SpMat tgt = ow ? src : as<MSpMat>( clone( m ) );
    F<double,ind,ind> func = get_func( static_cast<Weight>( mode ), N, rmrg, cmrg, positive );

    // CoocWorker wrkr( src, tgt, func );
    // RcppParallel::parallelFor( 0, src.outerSize(), wrkr );
    for( ind i = 0; i < src.rows(); ++i ) {
        for( SpInIt srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
            double v = srcIt.value();
            ind r    = srcIt.row();
            ind c    = srcIt.col();
            // tgtIt.valueRef() = func( v, r, c );
            tgtIt.valueRef() = -1;
        }
    }

    S4 out = wrap( tgt );
    out.slot( "Dimnames" ) = m.slot( "Dimnames" );
    return out;
}


