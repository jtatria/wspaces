#' Build a factor from max/min values in a series of variables
#'
#' @param d   A data frame
#' @param var An optional list of variables to compare
#' @param min If true, use minimum instead of maximum
#'
#' @return A factor of the same length as d, where the value for each observation will correspond
#'         to the name of the column in var that
#'         contains the minimum or maximum value for that row.
#'
#' @examples
#' d <- as.data.frame( matrix( rnorm( 300 ), nrow = 100, ncol = 3 ) )
#' names( d ) <- c( 'red', 'blue', 'green' )
#' count_to_factor( d )
#'
#' @export
#' @importFrom magrittr "%>%"
counts_to_factor <- function( d, var, min = FALSE, conf = FALSE ) {
  if( missing( var ) ) var = names( d )
  if( !all( apply( d[ , var ], 2, is.real ) ) ) stop( "Count variables must be numeric or integer" )
  return( factor( apply( d[,var], 1, findMax ), var ) )
}

#' Compute a (P)PMI matrix from raw frequency counts, stored as a column-oriented sparse matrix.
#'
#' @param m    A colun-oriented sparse matrix containing cooccurrence counts
#' @param ppmi A logical value indicating whether negative values should be truncated to 0 (i.e.
#'             compute PPMI instead of PMI). TRUE by default.
#' @param ow   A logical value indicating wether the result should be destructively copied over the
#'             input matrix. TRUE by default.
#'
#' @return The PMI value for each cell is equal to log( p(i,j) / p(i)p(j) ), i.e. the log of the
#'         observed probability over the expected probability. The PPMI truncates negative values
#'         to 0. WARNING: In order to maintain sparsity, -Inf values in the non-truncated case are
#'         replaced by 0.
#'
#' @export
cooc_to_pmi <- function( m, ppmi = TRUE, ow = TRUE ) {
  if( ncol( m ) + 1 != length( m@p ) ) {
    stop( "Wrong format for (P)PMI computation: m is not column-oriented compressed." )
  }

  # marginals
  N = sum( m )
  rp = Matrix::rowSums( m ) / N
  cp = Matrix::colSums( m ) / N

  # log function to use
  f <- ifelse( ppmi, log, log1p )

  # output vector
  out = vector( mode = 'numeric', length = length( m@x ) )

  # m@p contains offsets into m@i and m@x for each column
  for( col in 1:length( m@p ) ) {
    # i+1 because m@p always start with a 0. m@i[] + 1, because Matrix uses 0-indices internally.
    # row indices for this column
    rows = m@i[ col:m@p[ col + 1 ] ] + 1
    # values for this column
    vals = m@x[ col:m@p[ col + 1 ] ]
    # compute and assign for all rows in this column
    out[ col:m@p[ col + 1 ] ] <- f( ( vals / N ) / rp[rows] * cp[col] )
  }

  # overwrite m or return new
  if( ow ) {
    m@x <- out
    return( m )
  } else {
    return( Matrix::sparseMatrix( i = m@i, p = m@p, x = out ) )
  }
}

findMax <- function( x, ties = NA ) {
  i = which( x == max( x ) )
  i = ifelse( length( i ) > 1, ties, i )
  ifelse( is.null( names( x ) ), i, names( x )[i] )
}

is.real <- function( x ) {
  is.numeric( x ) || is.integer( x )
}

