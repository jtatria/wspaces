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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>

extern "C" {
  #include "SparseMatrix.h"
}

using namespace Rcpp;

//' @importClassesFrom Matrix dgTMatrix
S4 triplets( S4 spm ) {
    if( !spm.is( "Matrix" ) ) stop( "spm is not a sparseMatrix" );
    if( spm.is( "dgTMatrix" ) ) return spm;
    IntegerVector i = spm.slot("i");
    IntegerVector p = spm.slot("p");
    NumericVector x = spm.slot("x");
    Function f = Environment::namespace_env( "Matrix" )["sparseMatrix"];
    return f( _["i"]=i, _["p"]=p, _["x"]=x, _["index1"]=false, _["giveCSparse"]=false );
}

//' @importClassesFrom Matrix dgCMatrix
S4 compressed( S4 spm ) {
    if( !spm.is( "Matrix" ) ) stop( "spm is not a sparseMatrix" );
    if( spm.is( "dgCMatrix" ) ) return spm;
    IntegerVector i = spm.slot("i");
    IntegerVector j = spm.slot("j");
    NumericVector x = spm.slot("x");
    Function f = Environment::namespace_env( "Matrix" )["sparseMatrix"];
    return f( _["i"]=i, _["j"]=j, _["x"]=x, _["index1"]=false, _["giveCSparse"]=true );
}

//' @importClassesFrom Matrix dgTMatrix
S4 dgTMatrix( IntegerVector i, IntegerVector j, NumericVector x ) {
    Function f = Environment::namespace_env("Matrix")["sparseMatrix"];
    return f( _["i"]=i, _["j"]=j, _["x"]=x, _["index1"]=false );
}

//' @importClassesFrom Matrix dgCMatrix
S4 dgCMatrix( IntegerVector i, IntegerVector p, NumericVector x ) {
    Function f = Environment::namespace_env("Matrix")["sparseMatrix"];
    return f( _["i"]=i, _["p"]=p, _["x"]=x, _["index1"]=false );
}

//' Load sparse matrices from disk.
//'
//' Reads a sparse matrix from a binary file containing a sequence of (int,int,double) triplets.
//'
//' @param file A character vector indicating the path to the file to be read.
//'
//' @return A Matrix::sparseMatrix built from the triplets read from the given binary file.
//'
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
S4 load_spm( std::string file ) {
    SparseMatrix* spm = spm_read( file.c_str() );
    double* x = spm->entries;
    int*    i = spm->row_indices;
    int*    j = spm->col_indices;
    IntegerVector iv( i, i + spm->entry_count );
    IntegerVector jv( j, j + spm->entry_count );
    NumericVector xv( x, x + spm->entry_count );
    return dgTMatrix( iv, jv, xv );
}

//' Save sparse matrices to disk.
//'
//' Writes a sparse matrix object to a binary file as a sequence of (int,int,double) triplets.
//'
//' @param m    A sparseMatrix instance
//' @param file A file name
//'
// [[Rcpp::export]]
void save_spm( S4 m, std::string file ) {
    m = triplets( m );
    IntegerVector i( m.slot( "i" ) );
    IntegerVector j( m.slot( "j" ) );
    NumericVector x( m.slot( "x" ) );
    SparseMatrix* spm;
    spm->entry_count = x.size();
    spm->row_indices = i.begin();
    spm->col_indices = j.begin();
    spm_write( spm, file.c_str() );
    spm_free( spm );
}
