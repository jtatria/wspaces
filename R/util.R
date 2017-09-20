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

#' Kullback-Leibler divergence
#'
#' The Kullback-Leibler divergence of p from q is equal to
#' \eqn{ \sum_{i} P(i)\ log\ { \frac{P(i)}{Q(i)} } }.
#' KL is valid only when q_i = 0 iff p_i = 0. If p or q contain any zeros, a correction will be
#' applied, as per the 'zeros' parameter.
#'
#' If either p or q contain any zeros, the value of the 'zeros' parameter determines what action
#' to take, before computing x / sum( x ):
#' 'drop' will drop all with zero values in p or q:
#' \code{ x <- x[ p != 0 & q !=0 ] }.
#' 'smooth' will smooth both vectors with a Dirichlet prior:
#' \code{ x <- ( x + 1 ) / ( sum( x ) + length( x ) ) }.
#' 'convex' will apply a convex combination of x with a uniform distribution:
#' \code{ x <- x + ( 1 / length( x ) ) }.
#' 'bail' will error out.
#'
#' @param p     A frequency or probability vector. Coerced to prob by p / sum( p )
#' @param q     A frequency or probability vector. Coerced to prob by q / sum( q )
#' @param zeros Strategy for dealing with zeros in p or q. See details.
#' @param base Logarithm base. Uses \eqn{e} by default, yielding measure in nats. Use 2 for measure
#'             in bits/shannons.
#'
#' @return A scalar value equal to the Kullback-Leibler divergence of p from q.
#'
#' @export
div_kulback_leibler <- function(
    p, q, base=exp( 1 ), zeros=c('drop','smooth','convex','bail'), conexd=rexp
) {
    if( length( p ) != length( q ) ) stop( 'vector lengths differ' )
    zero = match.arg( zero )
    if( any( p == 0 ) || any( q == 0 ) ) {
        if( zeros != 'bail' ) warning( paste( 'applying correction for 0 values:', zeros ) )
        switch( zeros,
            drop   = {
                valid <- ( p != 0 & q != 0 )
                p <- p[valid]
                q <- q[valid]
            },
            # I'm almost positive that these are exactly the same...
            smooth = {
                p <- ( p + 1 ) / sum( p ) + length( p )
                q <- ( q + 1 ) / sum( q ) + length( q )
            },
            convex = {
                p <- p + ( 1 / length( p ) )
                q <- q + ( 1 / length( q ) )
            },
            bail = stop( '0 values in either p or q' )
        )
    }
    p <- p / sum( p )
    q <- q / sum( q )
    return( sum( p * log( p / q, base=base ) ) )
}

#' Jensen-Shannon divergence
#'
#' Computes the Jensen-Shannon smoothing reflection of the Kullback-Leibler divergence.
#'
#' The Jensen-Shannon divergence is equal to the arithmetic mean between the Kullback-Leibler
#' divergences of p from m and q from m, where m is equal to the elemen-wise mean between p and q.
#'
#' This function only computes the element-wise mean between p and q, and the mean between the KL
#' of p and from m. The actual Kullback-Leibler divergences are computed by
#' \code{div_kulback_leibler(p,q)}.
#'
#' Unlike the underlying Kullback-Leibler, this measure is reciprocal.
#'
#' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
#' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
#' @param ... further arguments passed to the underlying Kullback-Leibler implementation
#'            (@seealso div_kulback_leibler)
#' @param base Logarithm base. Uses \eqn{e} by default, yielding measure in nats. Use 2 for measure
#'             in bits/shannons.
#'
#' @return a scalar value equal to the Jensen-Shannon divergence between p and q.
#'
#' @export
div_jensen_shannon <- function( p, q, ... ) {
    m <- .5 * ( p + q )
    kl_p <- div_kulback_leibler( p, m, ... )
    kl_q <- div_kulback_leibler( q, m, ... )
    return( .5 * ( kl_p + kl_q ) )
}

#' Bhattacharyya divergence
#'
#' Computes the Bhattacharyya divergence between the probability distributions represented by the
#' given frequency or probabality vectors. The Bhattacharyya divergence is equal to the reciprocal
#' of the natural logarithm of the Bhattacharyya coefficient between p and q:
#' \eqn{ -log( \sum_{i}^{n} \sqrt{ p_i*q_i } ) }
#'
#' This function just computes the natural logarithm reciprocal. The coefficient is computed by
#' \code{bha_coef( p, q )}.
#'
#' This measure is not a proper distance because it does not respect the triangle inequality.
#' If a proper distance is needed, use the Hellinger distance (\code{dist_hellinger( p, q )}), also
#' based on the Bhattacharyya coefficient.
#'
#' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
#' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
#'
#' @return a scalar value equal to the Bhattacharyya divergence between p and q.
div_bhattacharyya <- function( p, q ) {
    return( -1 * log( coef_bhattacharyya( p, q ) ) )
}

#' Hellinger distance
#'
#' Computes the Hellinger distance between the probability distributions represented by the given
#' frequency or probability vectors p and q. The Hellinger distance is equal to the square root of
#' the 1-complement of the Bhattacharyya coefficient between p and q:
#' \eqn{ \sqrt{ 1 - \sum_{i}^{n} \sqrt{ p_i*q_i } } }.
#'
#' This function only computes the square root of the 1-complement. The coefficient is computed by
#' \code{coef_bhattacharyya(p, q)}.
#'
#' Unlike the closely related Bhattacharyya divergence, this measure is a proper distance as it
#' respects the triangle inequality.
#'
#' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
#' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
#'
#' @return a scalar value equal to the Hellinger divergence between p and q.
#'
#' @export
dist_hellinger <- function( p, q ) {
    return( sqrt( 1 - coef_bhattacharyya( p, q ) ) )
}

#' Bhattacharyya coefficient
#'
#' Computes the Bhattacharyya overlap coefficient between the probability distributions represented
#' by the given frequency or probability vectors. The Bhattacharyya coefficient is equal to
#' \eqn{\sum_{i}^{n} \sqrt{ p_i*q_i } }.
#'
#' The given vectors will be coerced to probabilities by \code{x = x / sum( x )}.
#'
#' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
#' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
#'
#' @return a scalar value equal to the Bhattacharyya overlap coefficient between p and q.
#'
#' @export
coef_bhattacharyya <- function( p, q ) {
    p <- p / sum( p )
    q <- q / sum( q )
    return( sum( sqrt( p * q ) ) )
}
