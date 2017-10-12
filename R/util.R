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


#' Real-valued vectors
#'
#' wrapper for ( is.numeric( x ) || is.integer( x ) )
#'
#' @param x a vector.
#'
#' @return TRUE if x can be interpreted as a real number, FALSE otherwise.
#'
#' @export
is.real <- function( x ) {
    is.numeric( x ) || is.integer( x )
}

#' Generate rank classes with equal population size.
#'
#' Generates classes of approximately equal size over intervals of the rankings in x. i.e. splits
#' the vector of ranks in x in classes s.t. each class contains a similar number of observations.
#'
#' This is similar to the "Jenks" classification strategy.
#'
#' @param x       A sortable vector.
#' @param k       The desired number of classes.
#' @param factor  If false, return a numeric vector instead of a factor.
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
    clz <- classify(
        1:length( x ),
        seq( 0, length( x ), length( x ) / k ),
        factor=factor, labels=labels
    )
    return( clz[ order( srt ) ] )
}

#' Generate classes with equal mass.
#'
#' Generates classes from a vector of masses (frequencies) s.t. that each generated class has
#' equivalent total mass (and most likely very different sizes for highly skewed distributions like
#' e.g. a lexicon's tf).
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
    clz <- classify(
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
#'
#' @return A vector of length equal to \code{length( x )} with the interval in k for each value of x
#'
classify <- function( x, k, factor=TRUE, labels=NA ) {
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

#' Ratio between sets.
#'
#' Computes a ratio of the given aggregation function over the elements of the given vector for
#' members of the given numerator set and the elements of the given vector for the given denominator
#' set.
#'
#' This function can be used to implement many different measures of similarity/divergence not
#' included in this package. It is also used for community/vertex weighting function in the wspaces'
#' graph module.
#'
#' @param x     A vector from where to select inputs to the given aggregation function.
#' @param nset  A vector representing a subset of x over which to apply the given function for the
#'              numerator value.
#' @param dset  A vector representing a subset of x over which to apply the given function for the
#'              denominator value.
#' @param nfunc An aggregation function like sum or length for the numerator set. Defaults to
#'              length.
#' @param dfunc An aggregation function like sum or length for the denominator set. Defaults to
#'              nfunc.
#'
#' @return A scalar value equal to the value of the given function applied to the elements of nset
#'         over the value of the given funciton applied to the elements of dset.
#'
#' @export
#' @importFrom magrittr %>%
div_setratio <- function( x, nset, dset, nfunc=length, dfunc=nfunc ) {
    if( is.logical( nset ) && length( nset ) != length( x ) )
        warning( "Recycling short logical vector for nset" )
    if( is.logical( dset ) && length( dset ) != length( x ) )
        warning( "Recycling short logical vector for dset" )
    if( is.integer( nset ) && max( nset ) > length( x ) )
        stop( 'Invalid index vector for nset (max nset oob)' )
    if( is.integer( dset ) && max( dset ) > length( x ) )
        stop( 'Invalid index vector for dset (max nset oob)' )
    return( ( x[ nset ] %>% nfunc ) / ( x[ dset ] %>% dfunc ) )
}

