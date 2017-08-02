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

// [[Rcpp::depends(Matrix)]]
#ifndef IOTOOLS_
#define IOTOOLS_ 1

#include <Rcpp.h>
extern "C" {
  #include "SparseMatrix.h"
}

Rcpp::S4 dgTMatrix( Rcpp::IntegerVector i, Rcpp::IntegerVector j, Rcpp::NumericVector x );
Rcpp::S4 dgCMatrix( Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::NumericVector x );

#endif // IOTOOLS_
