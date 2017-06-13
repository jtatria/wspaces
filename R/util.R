#' Build a factor from max/min values in a series of variables
#'
#' @param d A data frame
#' @param var An optional list of variables to compare
#' @param min If true, use minimum instead of maximum
#'
#' @return A factor of the same length as d, where the value for each observation will correspond to the name of the column in var that
#' contains the minimum or maximum value for that row.
#'
#' @examples
#' d <- as.data.frame( matrix( rnorm( 300 ), nrow = 100, ncol = 3 ) )
#' names( d ) <- c( 'red', 'blue', 'green' )
#' count_to_factor( d )
#'
#' @export
#' @importFrom magrittr "%>%"
count_to_factor <- function( d, var, min = FALSE, conf = FALSE ) {
  if( missing( var ) ) var = names( d )
  if( !all( apply( d[ , var ], 2, is.real ) ) ) stop( "Count variables must be numeric or integer" )
  return( factor( apply( d[,var], 1, findMax ), var ) )
  # out <- d[ , var ] %>% apply( 1, function( x ) { # apply( 1 ): over rows
  #   k = ifelse( min, var[ x == min( x ) ], var[ x == max( x ) ] ) # name of col vector with max/min count
  #   return( ifelse( length( k ) == 0, NA, k ) ) # I think we needed this to deal with collisions
  # } ) %>% unlist %>% as.factor() # we want a naked vector, as factor
  # return( out )
}

#' Compute a (P)PMI matrix from raw frequency counts.
#'
#' @param m A (typically sparse) matrix containing cooccurrence counts
#' @param positive A logical value indicating whether negative values should be truncated to 0 (i.e. compute PPMI instead of PMI)
#'
#' @return The PMI value for each cell is equal to log( p(i,j) / p(i)p(j) ), i.e. the log of the observed probability over the
#' expected probability. The PPMI truncates negative values to 0. WARNING: In order to maintain sparsity, -Inf values in the
#' non-truncated case are replaced by 0.
#'
#' @export
#' @importFrom Matrix "sparseMatrix"
#' @importMethodsFrom Matrix "rowSums" "colSums"
count_to_pmi <- function( m, positive = TRUE ) {
  x <- m@x       # data
  N <- sum( m )  # total observations

  obs <- x / N
  exp <- ( m != 0 ) * ( ( rowSums( m ) / N ) * colSums( m ) / N )

  if( positive ) {
    m@x <- log1p( obs / exp@x )
  } else {
    warning( "Using raw log instead of log1p. 0 entries in output should be -Inf but are set to 0 to maintain sparsity" )
    m@x <- log( obs / exp@x )
  }

  return( m )
}

findMax <- function( x, ties = NA ) {
  i = which( x == max( x ) )
  i = ifelse( length( i ) > 1, ties, i )
  ifelse( is.null( names( x ) ), i, names( x )[i] )
}

is.real <- function( x ) {
  is.numeric( x ) || is.integer( x )
}

