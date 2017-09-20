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
        Implementation of additive and difference-weigthed cooccurrence retrieval models.
***************************************************************************************************/

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends( RcppParallel )]]
#include "tools.h"
#include <functional>
#include <RcppParallel.h>

using namespace Rcpp;

enum Sim {
    ADDITIVE  = 0,
    DWEIGHTED = 1,
    COSINE    = 2,
};

struct SimWorker : public RcppParallel::Worker {
    const Mat src;
    const bool symm;
    const bool self;
    std::function<double(Vec,Vec,int,int)> func;

    Mat tgt;

    SimWorker(
        const Mat src, const Mat tgt, bool self, bool symm,
         std::function<double(Vec,Vec,int,int)> func
    )
        : src( src ), tgt( tgt ), func( func ), self( self ), symm( symm ) {}

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            for( std::size_t j = 0; j < ( symm ? i + ( self ? 1 : 0 ) : src.rows() ); j++ ) {
                //if( self && i == j ) continue; // redundant. i + self ? 1 : 0 should cover it.
                double v = func( src.row( i ), src.row( j ), i, j );
                tgt( i, j ) = v;
                if( symm ) tgt( j, i ) = v;
            }
        }
    }
};

std::function<double(Vec,Vec,int,int)> get_sim( Sim s ) {
    std::function<double(Vec,Vec,int,int)> func;
    switch( s ) {
        case ADDITIVE  : break;
        case DWEIGHTED : break;
        case COSINE    : break;
        default : stop( "Unknown similarity measure requested" );
    }
    return func;
}

bool symmetric( Sim s ) {
    switch( s ) {
        case ADDITIVE  : return false;
        case DWEIGHTED : return false;
        case COSINE    : return true;
        default : stop( "Unknown similarity measure requested" );
    }
}

inline std::function<double(Vec,Vec,int,int)> additive() {
    return []( Vec vi, Vec vj, int i, int j ) -> double {
        double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
        double d   = ( ( ii > 0 && ji > 0 ? ii + ji : 0 ) + ( ij > 0 && jj > 0 ? ij + jj : 0 ) );
        double num = (
            ( vi.array().sign() * vj.array().sign() ).array() * ( vi + vj ).array()
        ).sum() - d;
        double den = ( vi.sum() - ( ( ii > 0 ? ii : 0 ) + ( ii > 0 ? ii : 0 ) ) );
        return num / den;
    };
}

std::function<double(Vec,Vec,int,int)> dweighted() {
     return []( Vec vi, Vec vj, int i, int j ) -> double {
        double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
        double d   = ( std::min( ii, ij ) + std::min( ji, jj ) );
        double num = vi.cwiseMin( vj ).sum() - d;
        double den = vi.sum() - ( ii + ji );
        return num / den;
    };
}

std::function<double(Vec,Vec,int,int)> cosine() {
    return []( Vec vi, Vec vj, int i, int j ) -> double {
        double num = ( vi.dot( vj ) );
        double den = ( vi.norm() * vj.norm() );
        return num / den;
    };
}

// TODO : add wrapper.

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
RMatD sim_dweighted( RMatD m, bool self=false ) {
    Mat out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, self, false, []( Vec vi, Vec vj, int i, int j ) -> double {
        double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
        double d   = ( std::min( ii, ij ) + std::min( ji, jj ) );
        double num = vi.cwiseMin( vj ).sum() - d;
        double den = vi.sum() - ( ii + ji );
        return num / den;
    } );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return wrap( out );
}

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
RMatD sim_additive( RMatD m, bool self=false ) {
    Mat out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, self, false, []( Vec vi, Vec vj, int i, int j ) -> double {
        double ii  = vi.coeff( i ), ij = vj.coeff( i ), ji = vi.coeff( j ), jj = vj.coeff( j );
        double d   = ( ( ii > 0 && ji > 0 ? ii + ji : 0 ) + ( ij > 0 && jj > 0 ? ij + jj : 0 ) );
        double num = (
            ( vi.array().sign() * vj.array().sign() ).array() * ( vi + vj ).array()
        ).sum() - d;
        double den = ( vi.sum() - ( ( ii > 0 ? ii : 0 ) + ( ii > 0 ? ii : 0 ) ) );
        return num / den;
    } );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return wrap( out );
}

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
RMatD sim_cosine( RMatD m, bool self=false ) {
    Mat out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, self, true, []( Vec vi, Vec vj, int i, int j ) -> double {
        double num = ( vi.dot( vj ) );
        double den = ( vi.norm() * vj.norm() );
        return num / den;
    } );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return wrap( out );
}


