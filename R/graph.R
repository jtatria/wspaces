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
####                 Functions for creating and manipulating semantic networks                  ####
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
#' This function will call \code{igraph::components( g )} if no value is given for the cmps
#' parameter. \emph{This can be extremely slow}, but so far I have not been able to find or produce
#' a faster component determination strategy.
#' If you know of one, please contact me (e.g. if you know how to compute matrix kernels in less
#' than O(n^3))
#'
#' @param g    An igraph graph.
#' @param tol  An integer indicating the maximum size of detached components. Defaults to 1.
#' @param cmps An igraph components object. Will be computed over g if none given.
#'
#' @return \code{TRUE} if the given graph has no detached component larger than tol, \code{FALSE}
#'         otherwise.
#'
#' @export
#' @importFrom igraph components
graph_connected <- function( g, tol=1, cmps=components( g ) ) {
    if( cmps$no == 1 ) return( TRUE )
    if( max( cmps$csize[ -which.max( cmps$csize ) ] ) > tol ) {
        return( FALSE )
    } else {
        return( TRUE )
    }
}

#' Community and neighbourhood contributions
#'
#' This function computes the contributions of either a vertex's neighbourhood to its enclosing
#' community or a vertex's community to its sorrounding neighbourhood.
#'
#' Contributions are computed for each vertes as a ratio of a measure function of a numerator
#' set equal to a vertex's non-community crossing edges over a measure function over a denominator
#' set equal to a) all of its incident edges in the case of cluster to neighbourhood contribution
#' and b) all of its enclosing community in the case of neighbourhood to cluster contributions.
#'
#' The mode parameters controls which direction is computed, 'vc' computes negihborhood to cluster
#' contribution, 'cv' computes cluster to neighbourhood contributions.
#'
#' The attr parameter gives the names of an edge attribute to use as input for the measure function.
#' If none is given, edges are simply counted for both sets.
#'
#' The aggr parameter allows specification of a different aggregation function if edge attributes
#' are used.
#'
#' @param g    An igraph graph.
#' @param c    An igraph communities object.
#' @param mode Characer vector indicating the direction of the measure; one of 'vc' or 'cv'. See
#'             details.
#' @param vset An optional set of vertices to compute the measure for. Defaults to V( g ) (i.e. all
#'             vertices).
#' @param attr The name of an edge attribute to use as input for the aggregation funciton. Set to
#'             NA to use edge counts instead. Defaults to 'weight' if g is weighted.
#' @param aggr The aggregation function to use over both sets. Only applicable to edge attributes;
#'             ignored if attr is NA. Defaults to sum.
#'
#' @return A vector of length equal to the number of vertices in the given vset with the requested
#'         contribution measure values for each vertex in its assigned community.
#'
#' @export
#' @importFrom igraph membership as_adj
graph_cluster_contrib <- function(
    g, c, mode=c('cv','vc'), vset=V( g ),
    weight='weight', aggr=ifelse( is.na( attr ), length, sum )
) {
    chk_igraph( g ); chk_comm( c ); mode=match.arg( mode );
    k   <- igraph::membership( c )
    adj <- igraph::as_adj( g, type='both', attr=weight, sparse=FALSE )
    out <- switch( match.arg( mode ), cv=c2v_contrib( adj, k ), vc=v2c_contrib( adj, k ) )
    return( out )
}

#' Compute vertex weights for all clusters in the given communities object.
#'
#' This function is a simple wrapper over a vertex-cluster weight function s.t. the resulting
#' matrix will contain a a column for each cluster with a vector of vertex weights for all vertices
#' in the graph. All vertices will have a single non-zero entry in the column corresponding to its
#' assigned community.
#'
#' @param g            An igraph object.
#' @param cms          An igraph communities object.
#' @param contrib_func A function taking g and c as parameters to compute vertex weights. Defaults
#'                     to graph_cluster_contrib.
#' @param ...          Other parameters passed to the contrib_func.
#'
#' @return a matrix with as many rows as vertices in g and as many columns as cluster in cms.
#'
#' @export
#' @importFrom igraph membership
graph_cluster_contrib_matrix <- function( g, cms, contrib_func=graph_cluster_contrib, ... ) {
    chk_igraph( g ); chk_comm( cms )
    cwt  <- contrib_func( g, cms, ... )
    memb <- membership( cms )
    ks <- unique( memb )
    out <- matrix( NA, nrow=length( V( g ) ), ncol=length( ks ) )
    for( k in ks ) {
        out[,k] <- ifelse( memb == k, cwt, 0 )
    }
    return( out )
}

#' Cluster distances
#'
#' This function will compute a distance matrix for the clusters found in two different community
#' extraction results by first computing a vector of vertex weights for each cluster in each
#' solution and then computing a distance matrix between the clusters on both solutions based on
#' these vertex weight vectors.
#'
#' The contrib_func parameter must be a weighting function of the form f( g, c ) -> v, where g is
#' a graph, c is a communities object and the returned value v is a vector of weights for all
#' vertices in g. See graph_cluster_contrib for an example.
#'
#' The dist_func parameter must be a vector distance function of the form f( v1, v2 ) -> x, where
#' v1 and v2 are weight vectors, and x is a scalae value. See any of the similarity and divergence
#' functions in this package for examples.
#'
#' If the given cms1, cms2 objects come from different graphs (i.e. from different time periods,
#' etc), then both graphs must be supplied. If no second graph is supplied, it will be assumed
#' that both community solutions come from the same graph.
#'
#' The label_format parameter must be a format string containing two "\%d" format specs, no more
#' and no less. The first will be used for the cluster solution index (1 or 2), while the second
#' will be used for the cluster index within each solution. See the default value for an example.
#'
#' @param g1           An igraph object.
#' @param cms1         An igraph communities object, produced from g1.
#' @param cms2         An igraph communities object, produced from g1 or from g2 if given.
#' @param g2           An optional, second igraph object if cms2 is not from g1. NULL by default.
#' @param label_format A format sting to name rows and cols in the returned matrix. Defaults to
#'                     "c\%d_k\%d".
#' @param vid          Vertex attribute name for vertex ids, required if aligning across different
#'                     graphs.
#' @param fill         Logical indicating how to align cluster vertex vectors if aligning across
#'                     different graphs. TRUE (the default) means use vertex union. FALSE will use
#'                     vertex intersection.
#' @param contrib_func A function to compute vertex weights for a cluster.
#' @param dist_func    A function to compute distances between vertex weight vectors.
#' @param ...          Additional parameters passed to contrib_func.
#'
#' @return A matrix with as many rows as clusters in cms1 and as many columns as clusters in
#'         cms2 with the results of the pairwise application of dist_func.
#'
#' @export
graph_align_clusters <- function(
    g1, cms1, cms2, g2=NULL,
    label_format="c%d_k%d", vid="name", fill=TRUE,
    contrib_func=graph_cluster_contrib_matrix, dist_func=dist_hellinger, ...
) {
    filt <- if( !is.null( g2 ) ) {
        if( is.null( vertex_attr( g1, vid ) ) || is.null( vertex_attr( g2, vid ) ) ) {
            stop( "Invalid name for vertex id attr, needed for cross-graph alignment" )
        }
        if( fill ) { # union
            unique( c( vertex_attr( g1, vid ), vertex_attr( g2, vid ) ) )
        } else { # intersection
            vertex_attr( g1, vid ) %in% vertex_attr( g2, vid )
        }
    } else {
        V( g1  ) %>% as.integer()
    }

    g2 <- if( is.null( g2 ) ) g1 else g2
    m1 <- contrib_func( g1, cms1, ... )
    m2 <- contrib_func( g2, cms2, ... )
    rownames( m1 ) <- vertex_attr( g1, vid )
    rownames( m2 ) <- vertex_attr( g2, vid )

    out <- matrix( NA, nrow=ncol( m1 ), ncol=ncol( m2 ) )
    for( i in 1:nrow( out ) ) {
        for( j in 1:ncol( out ) ) {
            v1 <- m1[,i][filt] %>% ifelse( is.na( . ), 0, . )
            v1 <- m2[,j][filt] %>% ifelse( is.na( . ), 0, . )
            out[i,j] = dist_func( v1, v2 )
        }
    }

    if( !is.na( label_format ) ) {
        rownames( out ) <- sprintf( label_format, 1, 1:ncol( m1 ) )
        colnames( out ) <- sprintf( label_format, 2, 1:ncol( m2 ) )
    }

    return( out )
}

#' Create bipartite graph from a cluster alignment matrix.
#'
#' This function will transform a matrix containing the result of a cluster alignment call into a
#' bipartite graph.
#'
#' @param m A cluster alignment matrix.
#'
#' @return An igraph object containing a bipartite graph
#'
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix vertex_attr
graph_alignment_graph <- function( m ) {
    k1 = nrow( m ); k2 = ncol( m )
    gm <- matrix( 0.0, nrow=k1+k2, ncol=k1+k2 )
    rownames( gm ) <- c( rownames( m ), colnames( m ) )
    colnames( gm ) <- c( rownames( m ), colnames( m ) )
    gm[ 1:k1, (k1+1):ncol( gm ) ] <- m
    gm[ (k1+1):ncol( gm ), 1:k1 ] <- t( m )
    g <- graph_from_adjacency_matrix( gm, weighted=TRUE )
    vertex_attr( g, "type" ) <- c( rep( TRUE, k1 ), rep( FALSE, k2 ) )
    return( g )
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
    message( sprintf( sprintf( "[%%s]: |%%%ds|", w + 1 ),
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
