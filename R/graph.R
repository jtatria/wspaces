#' Prune a graph maintaining connectivity
#'
#' Removes edges from a graph such that the resulting graph maintains connectivity up to the given
#' tolerance.
#'
#' The given tolerance value indicates the maximum size of detached components that will be
#' tolerated when determining connectivity. i.e. a tolerance of one implies the graph will still be
#' considered connected if the size of any detached component after edge removal is no greater than
#' one.
#'
#' Edges are removed according to the significance order provided in the given edges vector.
#'
#' @param g An igraph graph.
#' @param edges A vector of edges, sorted according to their significance.
#' @param tol An integer indicating the maximum disconnected component allowed.
#' @param attr An edge attribute name, used to sort edges if no edge list provided. Ignored if
#'             edges is not NULL. Deafult 'weight' if it is.
#' @param desc Logical, indicating whether edges should be sorted in descending order according to
#'             attr. Ignored if edges is not NULL, default TRUE if it is.
#' @param dropV A logical indicating whether vertices in components below the given tolerance
#'              should be dropped from the graph. Default TRUE.
#' @return A subgraph of g with the least significant edges removed.
#'
#' @export
#'
#' @importFrom igraph V E edge_attr components
graph_prune_connected <- function( g, edges, tol = 1, dropV = TRUE, attr = 'weight', desc = TRUE ) {
    if( !inherits( g, 'igraph' ) ) stop( 'g is not a graph' )
    if( missing( edges ) ) {
      edges <- E( g )[ order( edge_attr( g, attr ), decreasing = desc ) ]
    }

    lo <- 1
    hi <- length( edges )
    while( hi - lo >= 1 ) {
        pivot <- ceiling( lo + ( hi - lo ) / 2 )
        if( graph_connected( g - edges[ pivot:length( edges ) ], tol = tol ) ) {
            if( hi == pivot ) lo <- pivot
            hi <- pivot
        } else {
            if( lo == pivot ) hi <- pivot
            lo <- pivot
        }
    }
    out <- g - edges[ hi:length( edges ) ]
    if( dropV ) {
      cmps <- components( out )
      out <- out - V( out )[ cmps$membership %in% which( cmps$csize <= tol ) ]
    }
    return( out )
}

#' Determine graph connectivity up to the given tolerance
#'
#' This function will return true if the maximum component size of all components detached from the
#' largest component is greater or equal to the given tolerance. E.g. a tolerance of one implies
#' that a graph will still be considered connected if the minor components are only orphan nodes.
#'
#' @param g An igraph graph.
#' @param tol An integer indicating the maximum size of detached components. Defaults to 1.
#'
#' @return \code{TRUE} if the given graph has no detached component larger than tol, \code{FALSE}
#'         otherwise.
#'
#' @export
#'
#' @importFrom igraph components
graph_connected <- function( g, tol = 1 ) {
    cmps <- components( g )
    if( cmps$no == 1 ) return( TRUE )
    if( max( cmps$csize[ -which.max( cmps$csize ) ] ) > tol ) {
        return( FALSE )
    } else {
        return( TRUE )
    }
}
