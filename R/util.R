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

####################################################################################################
#### Various utility functions for common data manipulation tasks.                              ####
####################################################################################################

# cooc_to_ppmi -------------------------------------------------------------------------------------
#' (Positive) Pointwise mutual information.
#'
#' Compute a (P)PMI matrix from raw frequency counts, stored as a column-oriented sparse matrix.
#'
#' The PMI value for each cell is equal to log( p(i,j) / p(i)p(j) ), i.e. the log of the observed
#' probability over the expected probability.
#'
#' The PPMI variant adds 1 to this value in order to truncate all values to 0 and maintain
#' sparsity.
#'
#' If rs or cs are given, the computation of the expected probabilities will use the values
#' contained in these vectors. This allows computing correct PMI values for partial cooccurrence
#' matrices, e.g. those built over a subsample of the entire corpus. If no marginal vectors are
#' given, these are computed from the row sums and column sums of the cooccurrence matrix,
#' respectively; this is equivalent to computing an 'empirical' PMI value from the given set of
#' observations, instead of computing it using the global corpus statistics.
#'
#' The computation of observed probabilities is obviously always based on the total number of
#' observations in the given cooccurrence matrix, independently of the marginal vectors used for
#' the expected probabilities.
#'
#' WARNING: In order to prevent a memory explosion, zero-values in the input matrix \emph{are never
#' calculated}, which has the numerical side-effect of replacing all -Inf values for the plain PMI
#' version with zeros. I have not yet thought of a better solution to distinguish -Infs from
#' actual 0s in the computed result.
#'
#' @param m     A colun-oriented sparse matrix containing cooccurrence counts.
#' @param rs    A numeric vector for row marginals (i.e. total DF for each term).
#' @param cs    A numeric vector for column marginals (i.e. total TF for each corpus segment).
#' @param ppmi  A logical value indicating whether negative values should be truncated to 0 (i.e.
#'              compute PPMI instead of PMI). TRUE by default.
#' @param ow    A logical value indicating wether the result should be destructively copied over
#'              the input matrix. FALSE by default.
#'
#' @return An (column-stored sparse) matrix, isomorphic to m_ with the (P)PMI values for m_.
#'         If ow == TRUE, m_ is replaced with this value.
#' @export
#' @importFrom Matrix rowSums colSums
# --------------------------------------------------------------------------------------------------
cooc_to_pmi <- function( m, rs = rowSums( m ), cs = colSums( m ), ppmi = TRUE, ow = FALSE ) {
    return( spm_pmi( m, rs, cs, ppmi, ow ) )
}

# make_codes ---------------------------------------------------------------------------------------
#' Build factor levels from combinations of variables.
#'
#' Constructs factor levels by assigning a value to combinations of columns. If mode is 'max', then
#' the columns will be treated as numeric and the associated factor level wiull be equal to the
#' name of the column contanining the maximum among the selected columns for each row. If mode is
#' 'comb', factor levels are constructed by assigning a character value to the concatenation of
#' values contained in the selected columns, optionally reduced to flags s.t. x = x > 0 if
#' flags==TRUE.
#'
#' @param d      A data frame.
#' @param ...    Arguments passed to \code{\link{dplyr::select}} to choose columns to compute
#'               levels on.
#' @param flags  Reduce values to booleans x = x > 0. Ignored if mode is not 'comb'. This will
#'               dramatically reduce the number of levels in the resulting factor. Default TRUE.
#' @param toint  Convert factors to integer values. Only makes sense when flags is TRUE, in which
#'               case the concatenation of logical values for x > 0 is interpreted as a binary
#'               integer.
#' @param target Character value indicating the name of the output factor. Default: 'code'
#' @param mode   One of 'max' or 'comb'. If 'max', value is equal to the row maximum for the
#'               selected columns. If 'comb', value is equal to the combination of all values in
#'               the selected column.
#' @param drop   Logical. If TRUE, \code{\link{droplevels}} is called on the factor before return.
#' @param ties   Function or NA. Function to break ties between max values. Ignored if mode is
#'               not 'max'.
#' @return a data frame with an extra variable containing a factor with the generated levels.
#'
#' @export
# --------------------------------------------------------------------------------------------------
make_codes <- function(
    d, ..., flags=TRUE, toint=TRUE, target='code', mode=c( 'max', 'comb' ), drop=TRUE, ties=NA
) {
    codify <- switch( match.arg( mode ),
        max  = function( x ) find_max( x, ties=ties ),
        comb = function( x ) ints_to_str( x, flags=flags, sep='' )
    )
    tmp <- dplyr::select_( d, .dots=lazyeval::lazy_dots( ... ) )
    if( !all( apply( tmp, 2, is.real ) ) & mode == max ) {
        stop( "Selected variables must be numeric or integer for max codes" )
    }
    d$codes_ <- apply( tmp, 1, codify( x ) )
    if( flags & toint ) d <- dplyr::mutate( d, codes_ = strtoi( codes_, base = 2 ) )
    lvls <- unique( d$codes_ )[ order( unique( d$codes_ ) ) ]
    d <- dplyr::mutate( d, codes_ = factor( codes_, levels=lvls ) )
    if( drop ) d <- dplyr::mutate( d, codes_ = droplevels( codes_ ) )
    dots_ <- list( 'codes_' ) %>% setattr( 'names', target )
    d %>% dplyr::rename_( .dots=dots_ )

    return( d )
}

# ints_to_str---------------------------------------------------------------------------------------
#' Paste a series of integer values into a character value.
#'
#' This function takes an integer vector as input and produces a string value equal to the
#' concatenation of the given vector, with the given separator added in between, and optionally
#' interpreting the given integer values as if they were boolean flags
#' \code{using as.integer( x > 0 )}.
#'
#' @param x     An integer vector
#' @param flags Logical. Convert values prior to concatenation by \code{as.integer( x > 0 )}.
#' @param sep   A character expresion tO use as separator between intger values. Empty by default.
#'
#' @return a character value containing the concatenation of x, optionally recoded as 0 and 1 and
#'         with the given separator in between values, if given.
#'
#' @export
# --------------------------------------------------------------------------------------------------
ints_to_str <- function( x, flags=TRUE, sep='' ) {
    if( flags ) x <- as.integer( x > 0 )
    out <- if( all( !is.na( x ) ) ) paste( x, collapse=sep ) else NA
    return( out )
}

# str_normalize ------------------------------------------------------------------------------------
#' Normalize strings.
#'
#' This function offers a consistent procedure for string normalization combining several common
#' operations into one call.
#'
#' @param x     A string (i.e. a character vector of length 1).
#' @param trim  Logical. Remove trailing and leading whitespace.
#' @param lower Logical. Reduce everything to lower case.
#' @param white Logical. Eliminate duplicate whitespace.
#' @param nl    Loigcal. Replace newlines with spaces.
#' @param punct Logicel. Remove punctuation.
#' @param word  Logical. Remove non-word characters.
#'
#' @return a reasonably normalized string.
#'
#' @examples
#' \code{
#' x <- "some Ugly    dirty\nstring!"
#' cat( x )
#' cat( normalizeString( x ) )
#' }
# --------------------------------------------------------------------------------------------------
str_normalize <- function( x, trim=TRUE, lower=TRUE, white=TRUE, nl=FALSE, punct = FALSE ) {
    if( !any( trim, lower, white, nl, punct ) ) return( x ) # All flags are false, nothing to do.
    if( lower ) x <- tolower( x )
    if( trim  ) x <- gsub( '^\\s+(.*?)\\s+$', '\\1', x, perl = TRUE )
    if( white ) x <- gsub( '\\s+', ' ', x )
    if( nl    ) x <- gsub( '\\n+', ' ', x )
    if( punct ) x <- gsub( '[[:punct:]]', '', x )
    if( word )  x <- gsub( '\\W+', '', x, perl = TRUE )
    return( x )
}

# random_spm ---------------------------------------------------------------------------------------
#' Create random sparse matrix.
#'
#' Creates a random sparse matrix with the given dimension, fill rate and value function.
#'
#' @param rows Integer vector of length 1. Number of rows.
#' @param cols Integer vector of length 1. Number of cols.
#' @param p    Numeric vector of length 1. 0 < p < 1; fill rate.
#' @param v    Numeric vector of length 1 or suplier function for values.
#'
#' @export
#' @importFrom Matrix Matrix
# --------------------------------------------------------------------------------------------------
random_spm <- function( rows = 100, cols = 100, p = 0.1, v = function() 1 ) {
    v <- if( !is.function( v ) ) function() v else v
    d <- sample( rbinom( rows*cols, 1, p ), rows*cols, TRUE )
    d <- d * v();
    m <- matrix( d, rows, cols, byrow = TRUE )
    return( Matrix( m, sparse=T ) )
}

find_max <- function( x, ties=NA ) {
    i = which( x == max( x ) )
    i = ifelse( length( i ) > 1, ties, i )
    ifelse( is.null( names( x ) ), i, names( x )[i] )
}

is.real <- function( x ) {
    is.numeric( x ) || is.integer( x )
}

#' @export
"%.%" <- function( ... ) paste( ..., sep ='' )
