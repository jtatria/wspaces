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
#### Functions for creating and manipulating semantic networks                                  ####
####################################################################################################

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
#' Edges are removed according to the significance order provided in the given edges vector or,
#' if none is provided, by the order induced by the given attribute.
#'
#' @param g     An igraph graph.
#' @param edges A vector of edges, sorted according to their significance.
#' @param tol   An integer indicating the maximum allowed disconnected component size.
#' @param attr  An edge attribute name, used to sort edges if no edge list provided. Ignored if
#'              edges is not NULL. Deafults to 'weight' if it is.
#' @param desc  Logical, indicating whether edges should be sorted in descending order according to
#'              attr. Ignored if edges is not NULL, defaults to TRUE if it is.
#' @param dropV A logical indicating whether vertices in disconnected components below the given
#'              tolerance tol should be dropped from the graph. Defaults to TRUE.
#'
#' @return A subgraph of g with the least significant edges removed.
#'
#' @export
#' @importFrom igraph V E edge_attr components
graph_prune_connected <- function(
    g, edges, tol=1, dropV=TRUE, attr='weight', desc=TRUE, verbose=FALSE
) {
    chk_igraph( g )
    if( missing( edges ) ) {
      edges <- E( g )[ order( edge_attr( g, attr ), decreasing = desc ) ]
    }
    head <- 1
    tail <- length( edges )
    while( tail - head >= 1 ) { # run until there are no more edges between tail and head
        pivot <- ceiling( head + ( tail - head ) / 2 ) # find middle edge between tail and head
        # test for connectedness after deleting all edges above pivot
        if( graph_connected( g - edges[ pivot:length( edges ) ], tol=tol ) ) {
            # if connected, set pivot to tail. if tail is aleady pivot, set head to pivot.
            if( verbose ) bins_prog( head, tail, pivot, length( edges ), lab='  connected  ' )
            if( tail == pivot ) head <- pivot
            tail <- pivot
        } else {
            # if not connected, set pivot to head. if head is aleady pivot, set tail to pivot.
            if( verbose ) bins_prog( head, tail, pivot, length( edges ), lab='not connected' )
            if( head == pivot ) tail <- pivot
            head <- pivot
        }
    }
    out <- g - edges[ tail:length( edges ) ]
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
#' that a graph will still be considered connected if the minor components are at most orphan nodes,
#' a tolerance of 2 implies the same if the minor components are at most only dyads, etc.
#'
#' @param g   An igraph graph.
#' @param tol An integer indicating the maximum size of detached components. Defaults to 1.
#'
#' @return \code{TRUE} if the given graph has no detached component larger than tol, \code{FALSE}
#'         otherwise.
#'
#' @export
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

#' Weight of ego's internal incident edges over all of its incident edges.
#'
#' Computes the ratio between ego's non-community crossing incident edges and all of its incident
#' edges. This measure is also known as ego's community's "contribution" to ego's neighbourhood.
#'
#' This can be interpreted as the 'weight' of ego's community on its inmediate neighbourhood.
#'
#' @param comm A communities object. The result of a cluster-finding function on g.
#' @param g The graph.
#' @param v A vector of vertices in g. Defualts to all vertices in g.
#' @param attr An edge attribute to sum over the two sets. Set to NA to use cardinalities instead.
#' @param aggr A function to combine values for the two sets. Defaults to sum.
comm_vertex_weight <- function( comm, g, v=V( g ), attr='weight', aggr=sum ) {
    chk_igraph( g ); chk_comm( comm )
    intra <- !crossing( comm, g )
    out <- vapply( v, function( v ) {
        # ego's intra-cluster incident edges
        nset <- ( ( E( g ) %in% incident( g, v ) ) & intra )
        # all of ego's incident edges
        dset <- ( ( E( g ) %in% incident( g, v ) ) )
        # edge attribute or edges themselves
        evec <- if( is.na( attr ) ) E( g ) else edge_attr( g, attr )
        # sum attribute or count edges.
        agrf <- if( is.na( attr ) ) length else sum
        return( set_ratio( evec, nset, dset, agrf ) )
    }, numeric( 1 ) )
    return( out )
}

#' Weight of ego's internal incident edges over all internal edges.
#'
#' Computes the ratio between ego's non-community crossing incident edges and all of its community's
#' internal edges. This measure is also known as ego's neighbourhood's "contribution" to its 
#' community.
#'
#' This can be interpreted as the 'weight' of ego's inmediate neighbourhood on its community.
#'
#' @param comm A communities object. The result of a cluster-finding function on g.
#' @param g The graph.
#' @param v A vector of vertices in g. Defualts to all vertices in g.
#' @param attr An edge attribute to sum over the two sets. Set to NA to use cardinalities instead.
#' @param aggr A function to combine values for the two sets. Defaults to sum.
vertex_comm_weight <- function( comm, g, v=V( g ), attr='weigth', aggr=sum ) {
    chk_igraph( g ); chk_comm( comm )
    intra <- !crossing( comm, g )
    memb <- membership( comm )
    out <- vapply( v, function( v ){
        # ego's intra-cluster incident edges
        nset <- ( ( E( g ) %in% incident( g, v ) ) & intra )
        # all intra-cluster edges for ego's cluster
        dset <- E( g )[ inc( V( g )[ memb == v ] ) ]
        # edge attribute or edges themselves
        evec <- if( is.na( attr ) ) E( g ) else edge_attr( g, attr )
        # sum attribute or count edges.
        agrf <- if( is.na( attr ) ) length else sum
        return( set_ratio( evec, nset, dset, agrf ) )
        return( num / den )
    }, numeric( 1 ) )
    return( out )
}

#' @importFrom igraph is.igraph
chk_igraph <- function( g ) {
    if( !igraph::is.igraph( g ) ) {
        stop( 'g is not an igraph graph!' )
    }
}

chk_comm <- function( comm ) {
    if( !'communities' %in% class( comm ) ) {
        stop( 'comm is not an igraph communities!' )
    }
}

bins_prog <- function( head, tail, pivot, range, lab, w=60 ) {
    lo <- ( ( ( 0     + head  ) / range ) * w ) %>% floor
    s1 <- ( ( ( pivot - head  ) / range ) * w ) %>% floor
    s2 <- ( ( ( tail  - pivot ) / range ) * w ) %>% floor
    hi <- ( ( ( range - tail  ) / range ) * w ) %>% floor
    message( sprintf( "[%s]: |%61s|",
        lab,
        paste( c(
            rep( ' ', lo ),
            rep( '-', s1 ),
            "|",
            rep( '-', s2 ),
            rep( ' ', hi )
        ), collapse='' )
    ) )
}
