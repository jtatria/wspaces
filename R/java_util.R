NS <- "edu.columbia.incite.obo"
OBO_CLZ <- paste( NS, 'OBOMain', sep = "/" )

#' Create a new OBO conf object
#'
#' The returned object may be used to modify corpus analysis parameters by using the native java
#' method 'set'. See examples.
#'
#' @return An OBOConf instance.
#'
#' @examples
#' conf <- make_obo_conf()
#' conf@set( 'wPre', 10 ) # Set cooccurrence window leading offset to 10.
#' conf$set( 'wPos', 10 ) # Set cooccurrence window trailing offset to 10.
#' obo <- make_obo( conf ) # obtain an OBO index interface object.
#'
#' @export
#' @importFrom rJava .jinit .jnew
make_obo_conf <- function() {
    .jinit()
    return( .jnew( paste( NS, 'OBOConf', sep = '/' ) ) )
}

#' Create OBO interface object.
#'
#' The returned OBO instance object will take its parameters from the given conf object. If no conf
#' is given, one will be created with default parameters.
#'
#' All index access functions require an OBO interface object as first parameter.
#'
#' @param conf An OBOConf instance, e.g. the value of \code{\link{make_obo_conf}}.
#'
#' @return An OBO index interface object.
#'
#' @export
#' @importFrom rJava .jinit .jnew %instanceof%
make_obo <- function( conf = make_obo_conf() ) {
    .jinit()
    if( !conf %instanceof% paste( NS, 'OBOConf', sep = "/" ) ) {
        stop( 'conf is not an OBOConf instance' )
    }
    return( .jnew( OBO_CLZ, conf ) )
}

#' Build document sets.
#'
#' Document sets define the sample of corpus segments that will be considered for all counting
#' functions, particularly \code{\link{get_ cooccurrences}}.
#'
#' If the given vector's length is greater than one, the filter will be constructed as the union of
#' all individual terms. If it is equal to one, it will be interpreted as a regular expression. If
#' its 0, the returned document set will have 0 documents.
#'
#' @param obo   An OBO interface object.
#' @param field A field over which to construct a document set.
#' @param terms A character vector containing terms to select documents for the DocSet.
#'
#' @return A DocSet instance that can be used to select observations from a corpus.
#'
#' @export
#' @importFrom rJava J .jinit .jarray %instanceof%
make_doc_set <- function( obo, field, terms ) {
    .jinit()
    if( !obo %instanceof% OBO_CLZ ) stop( 'Object is not an OBOMain instance' )
    if( length( field ) == 0 ) stop( 'empty field given' )
    if( length( field ) > 1 ) {
        warning( sprintf( 'multiple field names given for doc set. using %s', field[1] ) )
    }
    if( length( terms ) > 1 ) {
        return( obo$makeDocSet( field, .jarray( terms ) ) )
    } else {
        return( obo$makeDocSet( field, terms ) )
    }
}

#' Get all fields contained in the index.
#'
#' Obtains a character vector with all fields contained in the corpus index.
#'
#' The returned values may be used to e.g. list all terms in any field, construct document sets
#' over any one of the field's terms or combinations thereof, etc.
#'
#' @param obo   An OBO interface object.
#'
#' @return A character vector with the names of all fields present in the index.
#'
#' @export
#' @importFrom rJava .jinit J %instanceof%
get_index_fields <- function( obo ) {
    .jinit()
    if( !obo %instanceof% OBO_CLZ ) stop( 'Object is not an OBOMain instance' )
    return( J( paste( NS, 'util', 'RHelper', sep ='/'), 'fields', obo$indexReader() ) )
}

#' Get all terms in the given field.
#'
#' Obtains a character vector with all the terms found in the given field.
#'
#' Elements of the returned vector may be used to construct document sets over the same field
#' using \code{\link{make_doc_set}}.
#'
#' @param obo   An OBO interface object.
#' @param field A field name, e.g. an element of \code{\link{get_index_fields}}.
#'
#' @return A character vector containing all the terms found in the given field.
#'
#' @export
#' @importFrom rJava .jinit J %instanceof%
get_field_terms <- function( obo, field ) {
    .jinit()
    if( !obo %instanceof% OBO_CLZ ) stop( 'Object is not an OBOMain instance' )
    if( length( field ) > 1 ) {
        warning( sprintf( 'multiple field names given for doc set. using %s', field[1] ) )
    }
    return( J( paste( NS, 'util', 'RHelper', sep ='/'), 'terms', obo$indexReader(), field[1] ) )
}

#' Count cooccurrences over the given document set.
#'
#' Produces a sparse matrix containing cooccurrence counts for all terms in the analysis field over
#' all documents in the given document set.
#'
#' @param obo   An OBO interface object.
#' @param ds    A DocSet, built from \code{\link{make_doc_set}}.
#'
#' @return A sparse matrix containing cooccurrence counts for all terms in the given field in the
#' given sample of documents.
#'
#' @export
#' @importFrom Matrix sparseMatrix
#' @importFrom rJava .jinit
get_cooccurrences <- function( obo, ds = obo$defaultDocSet() ) {
    .jinit()
    cooc <- obo$countCooccurrences( ds )
    m <- cooc$arrays()
    return( sparseMatrix( i = m$i + 1, j = m$j + 1, x = m$x ) )
}
