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

//////////////////////////////// Similarity measures /////////////////////////////////////////

#include <stdio.h>

// [[Rcpp::depends( RcppParallel )]]
#include <RcppParallel.h>
#include "sim.hpp"

using namespace Rcpp;
// using namespace RcppParallel;

struct SimWorker : public RcppParallel::Worker {
    const Mat src;
    RcppParallel::RMatrix<double> tgt;
    const bool symm;
    const bool self;
    double (*func)(Vec,Vec,int,int);

    SimWorker(
        const Mat src, const NumericMatrix tgt, double(*func)(Vec,Vec,int,int), bool self, bool symm
    )
        : src( src ), tgt( tgt ), func( func ), self( self ), symm( symm ) {}

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            for( std::size_t j = 0; j < ( symm ? i : src.rows() ); j++ ) {
                if( self && i == j ) continue;
                double v = func( src.row( i ), src.row( j ), i, j );
                tgt( i, j ) = v;
                if( symm ) tgt( j, i ) = v;
            }
        }
    }
};

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
NumericMatrix sim_dweighted( NumericMatrix m, bool self=false ) {
    NumericMatrix out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, &simm_dw, self, false );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return out;
}

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
NumericMatrix sim_additive( NumericMatrix m, bool self=false ) {
    NumericMatrix out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, &simm_add, self, false );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return out;
}

//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
NumericMatrix sim_cosine( NumericMatrix m, bool self=false ) {
    NumericMatrix out( m.nrow(), m.ncol() );
    SimWorker wrkr( as<Mat>( m ), out, &simm_cos, self, true );
    RcppParallel::parallelFor( 0, m.nrow(), wrkr );
    return out;
}

