# IO functions to read and write wspaces data to/from disk

# entry sizes by difference to account for vector overhead.
comp_t <- object.size( c( complex( 1 ), complex( 1 ) ) ) - object.size( complex( 1 ) )
real_t <- object.size( c( 1.0, 1.0 ) ) - object.size( 1.0 )

#' @export
load_corpus <- function( dir, lexicon = "lexicon.tsv", frequencies = "frequencies.tsv", pos_table = "pos_counts.tsv", cooccur = "cooccur.bin" ) {
  lxcn = read_lexicon( lexicon )
  sfrq = read_frequencies( frequencies )
  posc = read_pos_counts( pos_table )
  cooc = read_cooccur( cooccur )

  return( list(
    "lexicon" = lxcn,
    "frequencies" = sfrq,
    "pos_counts" = posc,
    "cooccurrences" = cooc,
  ) )
}

#' @export
read_lexicon <- function( file, header = TRUE, sep = '@' ) {
  if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
  d <- read.table( file, quote = '', comment.char = '', sep = sep, header = header, row.names = 1 )
  return( d )
}

#' @export
read_frequencies <- function( file, header = TRUE, sep = '@', lxcn = NULL ) {
  if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
  d <- read.table( file, quote = '', comment.char = '', sep = sep, header = header, row.names = 1 )
  d[ is.na( d ) ] <- 0
  return( d )
}

#' @export
read_pos_counts <- function( file, header = TRUE, sep = '@' ) {
  if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
  d <- read.table( file, header = header, sep = sep, quote = "", comment.char = "", row.names = 1 )
  d[ is.na( d ) ] <- 0
  d$POS <- counts_to_factor( d )
  d$POS_conf <- apply( d[,-length( d )], 1, function( x ) max( x ) / sum( x ) )
  return( d )
}

#' @export
read_cooccur <- function( file ) {
  if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
  return( load_spm( file ) )
}

#' @export
#' @importFrom magrittr "%>%"
read_vectors <- function( filename, n = n, dim = dim, type = c( 'real', 'complex', 'binary' ), combine = FALSE ) {
  if( missing( type ) ) type = 'real'

  size <- as.integer( switch( type,
    real    = real_t,
    complex = comp_t,
    binary  = stop( "Binary vectors not supported" )
  ) )

  bytes <- file.info( filename )$size
  if( bytes %% ( 2 * size ) ) stop( "Truncated or corrupt file!" )

  if( !missing( n ) ) {
    dim <- ( bytes / size / 2 / n ) - 1
  } else if( !missing( dim ) ) {
    n <- ( bytes / size / 2 ) / ( dim + 1 )
  } else stop( "Either dim or n must be defined" )

  f <- file( filename, 'rb' )
  d <- f %>% readBin(
    switch( type, real = 'numeric', complex = 'complex', binary = stop( "Binary vectors not supported" ) ),
    n = ( dim + 1 ) * n * 2 ) %>% array( dim = c( dim + 1, n, 2 ) )
  close( f )

  if( combine ) {
    d <- d[,,1] + d[,,2]
  }

  return( d )
}

