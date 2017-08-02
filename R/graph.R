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

bins_prog <- function( head, tail, pivot, range, lab ) {
    w = 60
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

#' Scale edge attribute according to community structure.
#' 
#' Scales the given edge attribute by the given intra and inter community factors depending on 
#' whether the edges cross community boundaries or not.
#' 
#' @param edges  Edges in matrix form
#' @param mem    Vector of vertex membership
#' @param intraw Numeric factor for intra-community edges. Defaults to 1.
#' @param interw Numeric factor for inter-community edges. Defaults to 0.1 * intraw.
#' @param attr   Attribute scale. Defaults to 'weight'
#' 
#' @return A vector of scaled attributes attribute * intraw for edges in which the membership of 
#'         source and target is the same and attribute * interw if not.
#'         
#' @export
graph_comm_weight <- function( g, mem, intraw=1, interw=intraw*.1, attr='weight' ) {
    chk_igraph( g )
    w <- edge_attr( g, attr )
    edges <- get.edgelist( g )
    out <- vapply( 1:nrow( edges ), function( i ) {
        if( edge_crosses( edges[i,], mem ) ) {
            w[i] * interw
        } else {
            w[i] * intraw
        }
    }, 1.0 )
    return( out )
}

#' @importFrom igraph is.igraph
chk_igraph <- function( g ) {
    if( !igraph::is.igraph( g ) ) {
        stop( 'g is not an igraph graph!' )
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

cointet_loop <- function( m ) {
    out <- matrix( nrow=nrow( m ), ncol=ncol( m ) )
    for( i in 1:nrow( m ) ) {
        for( j  in 1:nrow( m ) ) {
            if( i != j ) {
                out[i,j] <- cointet_inner( m, i, j )
            }
        }
    }
    return( out )
}

cointet_inner <- function( m, w1, w2 ) {
    w1 <- if( is.character( w1 ) ) which( rownames( m ) == w1 ) else w1
    w2 <- if( is.character( w2 ) ) which( rownames( m ) == w2 ) else w2
    num <- sum( pmin( m[ w1, ][ -c( w1, w2 ) ], m[ w2, ][ -c( w1, w2 ) ] ) )
    den <- sum( m[ w1, ][ -c( w1, w2 ) ] )
    return( num / den )
}
