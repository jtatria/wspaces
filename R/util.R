#
# Copyright (C) 2017 José Tomás Atria <jtatria at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Utility functions for common data manipulation tasks.

#' Generate rank classes with equal size.
#'
#' Generates classes of approximately equal size over intervals of the rankings in x. i.e. splits
#' the vector of ranks in x in classes s.t. each class has similar size.
#'
#' @param x       A sortable vector.
#' @param k       The desired number of classes.
#' @param factor  If false, return a numerci vector instead of a factor.
#' @param labels  A vector of length k or a function to be applied over 1:k to be used as factor
#'                labels. Ignored if factor==FALSE
#' @param desc    Logical. Sort in descending order. TRUE by default.
#' @param no.sort Logical. Don't sort x, for cases in which x is already sorted (or results will
#'                be wrong).
#'
#' @return A vector of length equal to k with the class for each x.
#'
#' @export
class_rank <- function( x, k=100, factor=TRUE, labels=NA, desc=TRUE, no.sort=FALSE ) {
    x <- if( !no.sort ) x[ srt <- order( x, decreasing=desc ) ] else x
    clz <- mk_class(
        1:length( x ),
        seq( 0, length( x ), length( x ) / k ),
        factor=factor, labels=labels
    )
    return( clz[ order( srt ) ] )
}

#' Generate classes with equal mass.
#'
#' Generates classes from a vector of masses (frequencies) s.t. that each generated class has
#' equivalent mass (and most likely very different sizes).
#'
#' @param x      A vector of massese or frequencies.
#' @param k      The desired number of classes.
#' @param factor If false, return a numerci vector instead of a factor.
#' @param labels A vector of length k or a function to be applied over 1:k to be used as factor
#'               labels. Ignored if factor==FALSE.
#' @param log    Logical. Compute mass using log1p( x ) instead of x.
#' @param sort   Logical. Sort x before calculating mass and class; ensures members in each class
#'               have individual masses in the same order of magnitude.
#'
#' @return A vector of length equal to k with the class for each x.
#'
#' @export
class_mass <- function( x, k=100, factor=TRUE, labels=NA, log=FALSE, sort=TRUE, desc=TRUE ) {
    x <- if( sort ) x[ srt <- order( x, decreasing=desc ) ]
    x <- if( log ) cumsum( log1p( x ) ) else cumsum( x )
    clz <- mk_class(
        x,
        seq( 0, ( sup <- x[length( x )] ), sup / k ),
        factor=factor, labels=labels
    )
    return( if( exists( 'srt' ) ) clz[ order( srt ) ] else clz )
}

#' Wrapper for cut.
#'
#' Wrapper for cut( x, k ) that accepts vectors or a generator function for labels, and optionally
#' returns a numeric vector instead of a factor.
#'
#' @param x      Passed to cut as 'x'
#' @param k      Passed to cut as 'breaks'
#' @param factor Logical. Coerce to numeric if FALSE.
#' @param labels A vector of length equal to \code{length( k ) + 1} or a function to generate a
#'               vector of length \code{length( k ) + 1} from \code{1:lenght( k )} used as factor
#'               labels. Ignored if factor is FALSE.
#' @return A vector of length equal to \code{length( x )} with the interval in k for each value of x
#'
mk_class <- function( x, k, factor=TRUE, labels=NA ) {
    lbls <- NULL
    if( factor && !is.na( labels ) ) {
        if( is.vector( labels ) && ( length( labels ) == length( k ) - 1 ) ) lbls = labels
        else if( is.function( labels ) ) lbls = labels( 1:( k - 1 ) )
        else warning( 'Given labels is not function or has incorrect length.' )
    }
    clz <- cut( x, k, labels=lbls )
    if( !factor ) clz <- as.numeric( clz )
    return( clz )
}

#' Create random sparse matrix.
#'
#' Creates a random sparse matrix with the given dimensions, fill rate, value function and
#' (if given) respecting the given marginal distributions. Useful for bootstrapping tests for
#' sparse matrices.
#'
#' @param rows    Number of rows.
#' @param cols    Number of columns.
#' @param p       Numeric s.t. 0 < p < 1; fill rate (= probability of non-zero value)
#' @param v       Numeric vector of length 1 or suplier function for values. Ignored if row and col
#'                marginals are given.
#' @param row.mrg Numeric vector of length equal to rows with row marginals for computed matrices.
#'                If given, col.mrg must be given too.
#' @param col.mrg Numeric vector of length equal to cols with col marginals for computed matrices.
#'                If given, row.mrg must be given too.
#'
#' @return A Matrix::sparseMatrix.
#'
#' @export
#' @importFrom Matrix sparseMatrix
random_spm <- function( rows=100, cols=100, p=0.1, v=function() 1, row.mrg=NULL, col.mrg=NULL ) {
    # stop( 'not implemented yet' )
    # if( xor( is.null( row.mrg ), is.null( col.mrg ) ) ) {
    #     stop( 'Both margins needed when giving marginal values' )
    # }
    # v <- if( !is.null( row.mrg ) && !is.null( col.mrg ) ) {
    #     if( !is.vector( row.mrg, mode='numeric' ) || !is.vector( col.mrg, mode='numeric' ) ) {
    #         stop( 'marginals must be numeric vectors' )
    #     }
    #     if( length( row.mrg ) != rows || length( col.mrg ) != cols ) {
    #         stop( 'wrong dimensions for marginal vectors' )
    #     }
    #     function( i, j ) row.mrg[i] * col.mrg[j]
    # } else if( !is.function( v ) ) {
    #     function( i, j ) v
    # } else v
    v <- function( i, j ) 1
    nzs <- which( sample( rbinom( rows*cols, 1, p ), rows*cols, TRUE ) != 0 )
    i <- vector( 'integer', length=length( nzs ) )
    j <- vector( 'integer', length=length( nzs ) )
    x <- vector( 'numeric', length=length( nzs ) )
    for( nz in 1:length( nzs ) ) {
        i[nz] <- ( nz %/% rows ) + 1
        j[nz] <- ( nz %% rows ) + 1
        x[nz] <- v( i[nz], j[nz] )
    }
    browser()
    return( Matrix::sparseMatrix( i=i, j=j, x=x ) )
}

#' wrapper for ( is.numeric( x ) || is.integer( x ) )
is.real <- function( x ) {
    is.numeric( x ) || is.integer( x )
}
