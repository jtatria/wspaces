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

#' Jaccard coefficients and other set ratios.
#'
#' Computes a ratio of the given aggregation functions over the elements of the given vector for
#' members of the given numerator set and the elements of the given vector for members of the given
#' denominator set, i.e.: s1func( x[s1] ) / s2func( x[s2] ).
#'
#' This function can be used to implement many different measures of similarity/divergence not
#' included in this package. It is also used for cluster-vertex contribution scores in the wspaces'
#' graph module.
#'
#' If no second set is given, it is taken to be the universe, i.e. the entirety of x.
#'
#' If \code{jaccard} is TRUE, then the ratio is computed as per the Jaccard index:
#' \eqn{J = f(A \cap B) / f(A \cup B)}. If FALSE, the value is computed as a simple ratio
#' \eqn{f( A ) / f( B )}.
#'
#' @param x      A vector from where to select inputs to the given aggregation function.
#' @param s1     A logical or index vector representing a subset of x over which to apply the given
#'               function for the numerator value.
#' @param s2     A logical or index vector representing a subset of x over which to apply the given
#'               function for the denominator value. Defaults to \eqn{U}, i.e. all of x.
#' @param s1func An aggregation function like sum or length for the first set. Defaults to
#'               length.
#' @param s2fnuc An aggregation function like sum or length for the second set. Defaults to
#'               s1func, and there is seldom any reason to pass a different one.
#'
#' @param jaccard Logical. Compute a Jaccard distance. See details.
#'
#' @return A scalar value equal to the value of the given function applied to the elements of s1
#'         over the value of the given funciton applied to the elements of s2.
#'
#' @export
#' @importFrom magrittr %>%
div_set_ratio <- function( x, s1, s2=rep( TRUE, length( x ) ), s1func=length, s2func=s1func, jaccard=FALSE ) {
    if( is.logical( s1 ) && length( s1 ) != length( x ) ) {
        #warning( "Recycling short logical vector for s1" )
        stop( "logical vector is too short for s1" )
    }
    if( is.logical( s2 ) && length( s2 ) != length( x ) ) {
        #warning( "Recycling short logical vector for s2" )
        stop( "logical vector is too short for s2" )
    }
    if( is.integer( s1 ) && max( s1 ) > length( x ) ) {
        stop( 'Invalid index vector for s1 (max index oob)' )
    }
    if( is.integer( s2 ) && max( s2 ) > length( x ) ) {
        stop( 'Invalid index vector for s2 (max index oob)' )
    }
    if( jaccard ) {
        stop( 'Not implemented' ) # TODO
    }
    # TODO: hack pending https://github.com/igraph/igraph/issues/233
    s1 <- if( is.logical( s1 ) ) which( s1 ) else s1
    s2 <- if( is.logical( s2 ) ) which( s2 ) else s2
    # /hack
    return( ( s1func( x[ s1 ] ) ) / ( s2func( x[ s2 ] ) ) )
}

#' Extract marginals from global cooccurrence counts for sample cooccurrence counts
#'
#' Sparse matrices for corpus samples usually do not have the same dimension as global matrices,
#' because they do not store trailing zeros.
#'
#' Thus, the extraction of global marginals for e.g. weighting, requires resizing marginal vectors
#' to ensure that they are the same dimension as the sample matrix, i.e. as long as the extent of
#' the last row or column with non-zero values in the sample matrix, which will usually not be the
#' same as the last row or column in a global matrix.
#'
#' This function offers a measure of consistency for this operation for all corpus samples.
#'
#' @param s    A cooccurrence matrix for a corpus sample.
#' @param u    The global corpus cooccurrence matrix.
#' @param rank The rank to extract marginals on.
#' @param prob Logical. Transform marginal counts to a probability vector. Defults to TRUE.
#'
#' @return A row or column marginal vector with counts or probabilities of the correct dimension for
#'         use as a row or column marginal vector for weighting s
#'
#' @export
#' @importFrom Matrix rowSums colSums
sample_margins <- function( s, u, rank=c('row','col'), prob=TRUE ) {
    rank = match.arg( rank )
    mrg  <- if( rank == 'row' ) Matrix::rowSums( u ) else Matrix::colSums( u )
    fltr <- if( rank == 'row' ) rownames( s ) else colnames( s )
    mrg  <- mrg[ names( mrg ) %in% fltr ]
    mrg  <- if( prob ) mrg / sum( u ) else mrg
    return( mrg )
}

#' @rdname weight_cooc
#' @param u   The global corpus cooccurrence matrix.
#' @param ... Additional parameters passed to \link{weight_cooc}.
#' @export
weight_sample <- function( m, u=m, ... ) {
    r_mrg <- sample_margins( s, u, rank='row' )
    c_mrg <- sample_margins( s, u, rank='col' )
    wcooc <- weight_cooc( s, r_mrg, c_mrg, ... )
    return( wcooc )
}

