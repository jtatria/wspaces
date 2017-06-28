NS <- "edu.columbia.incite.obo"

#' @export
#' @importFrom rJava .jinit .jnew
make_obo_conf <- function() {
    .jinit()
    return( .jnew( paste( NS, 'OBOConf', sep = '.' ) ) )
}

#' @export
#' @importFrom rJava .jinit .jnew
make_obo <- function( conf ) {
    check_j( conf, paste( NS, 'OBOConf', sep = '.' ) )
    .jinit()
    return( .jnew( paste( NS, 'OBOMain', sep = '.' ), conf ) )
}

#' @export
#' @importFrom rJava J .jinit .jarray
make_doc_set <- function( obo, field, terms ) {
    check_j( obo, paste( NS, 'OBOMain', sep = '.' ) )
    .jinit()
    return( obo$makeDocSet( field, .jarray( terms ) ) )
}

#' @export
#' @importFrom rJava .jinit J
get_doc_fields <- function( obo ) {
    check_j( obo, paste( NS, 'OBOMain', sep = '.' ) )
    .jinit()
    clz <- J( paste( NS, 'uima', 'index', 'OBOLuceneFB', sep = '.' ) )
    return( clz )
}

#' @export
#' @importFrom Matrix sparseMatrix
get_cooccurrences <- function( obo, ds ) {
    check_j( obo, paste( NS, 'OBOMain', sep = '.' ) )
    cooc <- obo$countCooccurrences( ds )
    m <- cooc$arrays()
    return( sparseMatrix( i = m$i + 1, j = m$j + 1, x = m$x ) )
}

check_j <- function( obj, class ) {
    if( class( obj ) != 'jobjRef' ) {
        stop( 'object is not a java object.' )
    }
    clz <- gsub( "\\.", "/", class )
    if( obj@jclass != clz ) {
        stop( sprintf( "object is not an instance of %s", class ) )
    }
}
