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

#' @importFrom Matrix Matrix
#' @export
makeSpm <- function( rows = 100, cols = 100, p = 0.1, v = 1 ) {
  d <- sample( rbinom( rows*cols, 1, p ), rows*cols, TRUE )
  d <- d * v;
  m <- matrix( d, rows, cols, byrow = TRUE )
  return( Matrix( m, sparse=T ) )
}

findMax <- function( x, ties = NA ) {
  i = which( x == max( x ) )
  i = ifelse( length( i ) > 1, ties, i )
  ifelse( is.null( names( x ) ), i, names( x )[i] )
}

is.real <- function( x ) {
  is.numeric( x ) || is.integer( x )
}


