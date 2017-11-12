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

# RJava functions for creating and controlling OBO session objects.
# TODO: refactor to encapsulate native jav amethods in obo object.
# TODO: refactor to do without obo objects

NS       <- "edu/columbia/incite/obo"
OBO_CLZ  <- paste( NS, 'OBO', sep="/" )
CONF_CLZ <- paste( NS, 'OBOConf', sep='/' )
DS_CLZ   <- paste( NS, 'util', 'DocSet', sep='/' )
RHLP_CLZ <- paste( NS, 'util', 'RHelper', sep ='/' )

# TODO: add logic to check for destruction of JVM on rm of last obo object.

infof <- function( format, ... ) {
    .jinit
    msg <- sprintf( format, ... )
    log <- J( "java.util.logging.Logger", "getLogger", "edu.columbia.incite" )
    log$info( msg )
}

warnf <- function( format, ... ) {
    .jinit
    msg <- sprintf( format, ... )
    log <- J( "java.util.logging.Logger", "getLogger", "edu.columbia.incite" )
    log$warning( msg )
}

#' Create a new OBO conf object
#'
#' The returned object may be used to modify corpus analysis parameters by using the native java
#' method 'set'. See examples.
#'
#' @return An OBOConf instance.
#'
#' @examples
#' conf <- obo_mkconf()
#' conf@set( 'wPre', 10 ) # Set cooccurrence window leading offset to 10.
#' conf$set( 'wPos', 10 ) # Set cooccurrence window trailing offset to 10.
#' obo <- obo_new( conf ) # obtain an OBO index interface object.
#'
#' @export
#' @importFrom rJava .jinit .jnew
obo_mkconf <- function( ... ) {
    # TODO add par-like behaviour to conf
    args <- list( ... )
    .jinit()
    conf <- .jnew( CONF_CLZ )
    for( arg in names( args ) ) {
        par <- if( grep( "_(file|dir)", arg ) ) path.expand( args[[arg]] ) else args[[arg]]
        conf$set( arg, par )
    }
    return( conf )
}

#' Create OBO interface object.
#'
#' The returned OBO instance object will take its parameters from the given conf object. If no conf
#' is given, one will be created with default parameters.
#'
#' All index access functions require an OBO interface object as first parameter.
#'
#' @param conf An OBOConf instance, e.g. the value of \code{\link{obo_mkconf}}.
#'
#' @return An OBO index interface object.
#'
#' @export
#' @importFrom rJava .jinit .jnew
obo_new <- function( ..., conf=obo_mkconf( ... ) ) {
    .jinit()
    chk_clz( conf, CONF_CLZ )
    return( .jnew( OBO_CLZ, conf ) )
}

#' Rebuild corpus datasets.
#'
#' This function will rebuild all corpus datasets using the parameters currently found in the given
#' obo object's configuration.
#'
#' The produced data sets will be dumped into the current configuration's data directory, using the
#' configured filenames.
#'
#' @param obo An OBO interface object.
#' @param realod Logical indicating whether the new corpus files should be reloaded.
#' @param ...   Further arguments passed to load_corpus. Ignored if reload is FALSE.
#'
#' @return If reload, the value of load_corpus. NULL otherwise.
#'
#' @export
#' @importFrom rJava .jinit J
obo_rebuild_corpus <- function( obo, reload=TRUE, ... ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    J( obo, 'rebuildCorpus' )
    if( reload ) {
        load_corpus( obo$conf$dataDir()$toString(), ... )
    }
}

#' Build document sets.
#'
#' Document sets define the sample of corpus segments that will be considered for all counting
#' functions, particularly \code{\link{get_ cooccurrences}}.
#'
#' If the given vector's length is greater than one, the filter will be constructed as the union of
#' all individual terms. If it is equal to one, it will be interpreted as a regular expression. If
#' it's 0, the returned document set will be the empty set.
#'
#' @param obo   An OBO interface object.
#' @param field A field over which to construct a document set.
#' @param terms A character vector containing terms to select documents for the DocSet.
#'
#' @return A DocSet instance that can be used to select observations from a corpus.
#'
#' @export
#' @importFrom rJava .jinit J .jarray
obo_mkdocset <- function( obo, field=NULL, terms ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    # TODO: add support for int vectors?
    terms <- as.character( terms )
    field <- if( is.null( field ) ) obo$conf()$fieldDocId() else field
    if( length( field ) == 0 ) stop( 'empty field given' )
    if( length( field ) > 1 ) {
        infof( 'multiple field names given for doc set. using %s', field[1] )
    }
    if( length( terms ) > 1 ) {
        return( obo$makeDocSet( field, .jarray( terms ) ) )
    } else {
        return( obo$makeDocSet( field, terms ) )
    }
}

#' Intersect document sets.
#'
#' Combines the given documents sets such that the resulting document set is equal to their
#' intersection. This is equivalent to a logical AND over both document sets' vectors, but avoids
#' instantiating the full int vector.

#' @param ds1  A DocSet, built from \code{\link{obo_mkdocset}}.
#' @param ds2  A DocSet, built from \code{\link{obo_mkdocset}}.
#' @param free Logical indicating if the input sets shuld be destroyed in order to attempt recovery
#'             of their heap memory. Defaults to FALSE.
#'
#' @return A DocSet instance that can be used to select observations from a corpus equal to the
#'         intersection of the given doc sets.
#'
#' @export
#' @importFrom rJava .jinit J
obo_intersect <- function( ds1, ds2, free=FALSE ) {
    .jinit()
    chk_clz( ds1, DS_CLZ ); chk_clz( ds1, DS_CLZ )
    infof( "Intersecting doc set of size %d with doc set of size %d", ds1$size(), ds2$size() )
    ds <- ds1$intersect( ds2 )
    infof( "Resulting doc set has %d elements", ds$size() )
    if( free ) {
        J( RHLP_CLZ, "release", ds1 )
        J( RHLP_CLZ, "release", ds2 )
    }
    return( ds )
}

#' Obtain document sets for different document samples.
#'
#' Produces a doc set instance representing the requested sample. Samples are hardcoded for now,
#' pending a stale API to pass index queries back to the backend.
#'
#' Backend currently offers three samples: 'testimony' corresponds to all documents in a trial that
#' do not contain any legal entities. 'legal' corresponds to all documents in a trial that do
#' contain legal entities. 'trials' contain all documents in trial accounts.
#' \eqn{ testimony \cup legal = trials}.
#'
#' @param obo    An OBO interface object
#' @param sample The requested sample. One of "testimony", "legal", or "trials".
#'
#' @return A DocSet instance that can be used to select observations from a corpus corresponding
#'         to the requested sample.
#' @export
obo_sample <- function( obo, sample=c('testimony','legal','trials'), negate=FALSE ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    sample = match.arg( sample )
    ds <- if( negate ) {
        ds <- obo$complement( sample )
    } else {
        ds <- obo$getSample( sample )
    }
    return( ds )
}

#' Count cooccurrences over the given document set.
#'
#' Produces a sparse matrix containing cooccurrence counts for all terms in the analysis field over
#' all documents in the given document set.
#'
#' @param obo   An OBO interface object.
#' @param ds    A DocSet, built from \code{\link{obo_mkdocset}}.
#'
#' @return A sparse matrix containing cooccurrence counts for all terms in the given field in the
#' given sample of documents.
#'
#' @export
#' @importFrom Matrix sparseMatrix rowSums colSums
#' @importFrom rJava .jinit
obo_count_cooc <- function( obo, ds=obo$docSample(), shrink=TRUE ) {
    .jinit()
    chk_clz( obo, OBO_CLZ ); chk_clz( ds, DS_CLZ )
    prog <- J( RHLP_CLZ, "report" )
    cooc <- obo$countCooccurrences( ds, prog )
    ar <- cooc$arrays()
    m <- sparseMatrix( i=( ar$i + 1 ), j=( ar$j + 1 ), x=ar$x )
    J( RHLP_CLZ, "release", cooc )
    J( RHLP_CLZ, "release", ar )
    terms <- obo$lexicon()$arrays()$terms
    rownames( m ) <- terms[ 1:nrow( m ) ]
    colnames( m ) <- terms[ 1:ncol( m ) ]
    if( shrink ) m <- m[ rowSums( m ) > 0, colSums( m ) > 0 ]
    return( m )
}

#' Count frequencies over the given document set
#'
#' Produces a frequency table containing term frequencies for all terms in the lexicon in the
#' documents contained in the given doc set, split by the terms found in the configured split field.
#'
#' @param obo An OBO interface object
#' @param ds  A DocSet, built from \code{\link{obo_mkdocset}}. If none given, defaults to the given
#'            index's default doc sample.
#'
#' @return A lexical dataset containing term frequencies for the given documents over the
#'         configured field.
#'
#' @export
#' @importFrom rJava .jinit J
obo_count_freqs <- function( obo, ds=obo$docSample() ) {
    .jinit()
    chk_clz( obo, OBO_CLZ ); chk_clz( ds, DS_CLZ )
    prog <- J( RHLP_CLZ, "report" )
    off <- obo$conf()$freqFile()$toString()
    obo$conf()$set( obo$conf()$PARAM_FREQ_FILE, ( tmp <- tempfile() ) )
    freq <- obo$countFrequencies( ds, prog )
    obo$dumpFrequencies( freq )
    obo$conf()$set( obo$conf()$PARAM_FREQ_FILE, off )
    J( RHLP_CLZ, "release", freq )
    freqs <- read_frequencies( tmp )
    return( freqs )
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
#' @importFrom rJava .jinit J
obo_get_fields <- function( obo ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    return( J( RHLP_CLZ, 'fields', obo$indexReader() ) )
}

#' Get all terms in the given field.
#'
#' Obtains a character vector with all the terms found in the given field.
#'
#' Elements of the returned vector may be used to construct document sets over the same field
#' using \code{\link{obo_mkdocset}}.
#'
#' @param obo   An OBO interface object.
#' @param field A field name, e.g. an element of \code{\link{obo_get_fields}}.
#'
#' @return A character vector containing all the terms found in the given field.
#'
#' @export
#' @importFrom rJava .jinit J
obo_get_terms <- function( obo, field ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    if( length( field ) > 1 ) {
        warnf( 'multiple field names given for doc set. using %s', field[1] )
    }
    return( J( RHLP_CLZ, 'terms', obo$indexReader(), field[1] ) )
}

#' Get the total number of documents in the given index.
#'
#' @param obo An OBO interface object.
#'
#' @return An integer equal to the total number of documents in the given index.
#'
#' @export
#' @importFrom rJava .jinit
obo_numdocs <- function( obo ) {
    .jinit()
    chk_clz( obo, OBO_CLZ )
    return( obo$indexReader()$numDocs() )
}

#' @importFrom rJava %instanceof%
chk_clz <- function( obj, clz ) {
    if( !obj %instanceof% clz ) {
        stop( 'Object is not an instance of ' %.% clz )
    }
}

#' Get a copy of the lexicon.
#' @param obo An OBO interface object.
#' @return A lexicon data frame, with terms as row names and tf and df as columns.
#' @export
#' @importFrom rJava .jinit
obo_lexicon <- function( obo ) {
    .jinit()
    ar = obo$lexicon()$arrays()
    lxcn <- data.frame( row.names=ar$terms )
    lxcn[['tf']] <- ar$tf
    lxcn[['df']] <- ar$df
    return( lxcn )
}

#' Get all defined POS classes in the OBO corpus.
#' @return Character vector containing all POS class names in OBO.
#' @export
#' @importFrom rJava .jinit J
obo_all_pos <- function() {
    .jinit()
    return( J( RHLP_CLZ, 'posClasses' ) )
}

#' Get lexical POS classes in the OBO corpus.
#' @return Character vector containing lexical POS class names in OBO.
#' @export
#' @importFrom rJava .jinit J
obo_lex_pos <- function() {
    .jinit()
    return( J( RHLP_CLZ, 'lexClasses' ) )
}

