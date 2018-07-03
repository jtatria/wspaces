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

graph_make <- function(
    X, ls=rep( TRUE, nrow( X ) ), fs=rep( TRUE, ncol( X ) ), sim.func=NULL, 
    edge.mode=c( 'auto', 'directed', 'undirected' ), edge.normalize=TRUE, allow.loops=TRUE,
    prune.tol=1, prune.sort=NULL,
    vertex.data=NULL,
    cluster=TRUE, cluster.func=function( G ) igraph::cluster_louvain( igraph::as.undirected( G ) ),
    cluster.contribs=FALSE,
    ...
) {
    X %<>% X[ls,fs]
    if( !is.null( sim.func ) ) X %<>% sim.func
    edge.mode <- match.arg( edge.mode )
    edge.mode <- ifelse( edge.mode == 'auto',
        ifelse( isSymmetric( unname( X ) ),
            'undirected',
            'directed'
        ),
        edge.mode
    )
    
    G <- X %>% igraph::graph_from_adjacency_matrix(
        mode=edge.mode, diag=allow.loops, weighted=TRUE 
    )
    
    edges <- if( is.null( prune.sort ) ) {
        igraph::E( g )[ order( igraph::edge_attr( g, attr ), decreasing=TRUE ) ]
    } else if( is.function( prune.sort ) ) {
        igraph::E( g )[ order( prune.sort( g, ... ), decreasing=TRUE ) ]
    }
    G %<>% graph_prune( edges=edges, tol=prune.tol, desc=desc, dropV=TRUE, verbose=FALSE )
    
    if( edge.normalize ) {
        if( !is.null( igraph::edge_attr( g, 'weight' ) ) ) {
            w <- igraph::edge_attr( g, 'weight' )
            igraph::edge_attr( g, 'weight' ) <- w / max( w )
        } else {
            warning( 'Can\'t normalize weights in unweighted graph!' )
        }
    }
    
    if( !is.null( vertex.data ) ) G %<>% graph_add_vertex_data( vertex.data )
    
    if( cluster ) {
        K <- graph_cluster( G, cluster.func=cluster.func )
        G %<>% graph_add_cluster_data( G, K, contribs=cluster.contribs, quiet=FALSE )
    }
    
    obj <- structure( list(
        X=X, G=G, K=K,
        params=list(
            prune.sort     = prune.sort,
            sim.func       = sim.func,
            prune.tol      = prune.tol,
            edge.normalize = edge.normalize,
            allow.loops    = allow.loops
        )
    ), class='semnet' )
    return( g )
}

# Construction and pruning --------------------------------------------------------------------

#' #' Semantic Networks
#' #'
#' #' This function wraps a number of lower-level functions to generate semantic networks consistently
#' #' with a given set of parameters.
#' #'
#' #' In the simplest case, it will construct a first-order network using the given value of \code{X}
#' #' as an adjacency matrix, for all terms indicated by the given lexical sampling vector \code{ls}.
#' #'
#' #' If \code{x2_mode} is not \code{NA} or a value is given for \code{X2}, it will also compute a
#' #' higher order semantic network, either using the given value of \code{x2_mode} to compute a
#' #' similarity matrix calling \link{simdiv} or using the given \code{X2} value directly. If the
#' #' higher order matrix is computed internally, an optional \code{fs} vector can be used to select
#' #' a subset of available features, greatly speeding up computation at the cost of dropping some
#' #' features from consideration.
#' #'
#' #' In all cases this function will apply the relevant filtering for e.g. lexical sampling and
#' #' call \link{stitch} to construct an igraph graph object, then prune the resulting network up to
#' #' the given \code{tol} value by calling \link{prune} and detect communities over the pruned
#' #' network by calling the given \code{clust_func} clustering function.
#' #'
#' #' @param X       A cooccurrence matrix.
#' #' @param ls      A logical vector for lexical sampling. See \link{lexical_sample}.
#' #' @param fs      A logical vector for feature sampling. Defaults to all features.
#' #' @param tol     Integer. Connectivity tolerance for prunning. See \link{prune}.
#' #' @param x2_mode Integer. Mode for second-order network. See \link{simdiv} for possible values.
#' #'                Set to \code{NA} to skip higher-order network construction.
#' #' @param td      Optional term data to add as vertex attributes.
#' #' @param X2      A higher-order similarity matrix. Defaults to \code{NULL}. If present, a
#' #'                higher-order network will be constructed with this data and x2_mode will be
#' #'                ignored \emph{even if its set to \code{NA}}.
#' #' @param clust_function A community detection function. Defaults to \link{igraph::cluster_louvain}.
#' #' @param contribs Logical. Compute vertex-cluster contribution scores. Defaults to \code{FALSE}.
#' #'
#' #' @return If a higher order network is computed, a list containing two elements consisting of
#' #' objects of class "semnet" corresponding to the first- and second- order networks, named as
#' #' "order1" and "order2", respectively. If no higher-order network is requested
#' #' (i.e. \code{is.na( x2_mode ) && is.null( X2 ) } ), a single "semnet" class object containing the
#' #' network corresponding the \code{X} (i.e. the "order1" element in the list case).
#' #'
#' #' @export
#' graph_make <- function(
#'     X, ls, fs=rep( TRUE, ncol( X ) ), x2_mode=1, X2=NULL,
#'     tol=1, cluster_func=igraph::cluster_louvain, contribs=FALSE, td=NULL,
#'     quiet=FALSE
#' ) {
#'     S1 <- semnet_neu(
#'         X[ls,ls], tol=tol, cluster_func=cluster_func, contribs=contribs, td=td, quiet=quiet
#'     )
#' 
#'     if( !is.na( x2_mode ) || !is.null( X2 ) ) {
#'         if( !quiet ) message( 'Producing higher-order network; value will be list of length 2' )
#'         if( is.null( X2 ) ) {
#'             if( !quiet ) message( sprintf(
#'                 "Computing second order cooccurrences using mode %d. This may take a while...",
#'                 x2_mode
#'             ) )
#'             if( !quiet ) message( sprintf( 'Using feature vectors of length %d', sum( fs ) ) )
#'             X2 <- X[ls,fs] %>% simdiv( mode=x2_mode )
#'         }
#'         S2 <- X2 %>% semnet_neu(
#'             tol=tol, cluster_func=cluster_func, contribs=contribs, td=td, quiet=quiet
#'         )
#'     }
#' 
#'     if( is.null( S2 ) ) return( S1 )
#'     else return( list( order1=S1, order2=S2 ) )
#' }

semnet_neu <- function(
    X,tol=1, cluster_func=igraph::cluster_louvain, contribs=FALSE, td=NULL, quiet=FALSE
) {
    if( !quiet ) message( sprintf( "Stitching graph for %dx%d matrix with %f total tf mass",
        nrow( X ), ncol( X ), sum( X )
    ) )
    G <- X %>% graph_stitch( mode='undirected', vdf=td, tol=tol, verbose=!quiet )
    if( !quiet ) message( sprintf( "Clustering graph with %d vertices and %d edges",
        length( igraph::V( G ) ), length( igraph::E( G ) )
    ) )
    K <- G %>% graph_cluster( cluster_func=cluster_func, undirected=TRUE )
    G %<>% graph_add_cluster_data( K, contribs=contribs )
    obj <- list( X=X, G=G, K=K )
    class( obj ) <- 'semnet'
    return( obj )
}

#' Stitch network from the given adjacency matrix.
#'
#' This function will produce a full graph from all entries in the given adjacency matrix and then
#' prune the full graph by removing least significant edges while maintaining connectivity.
#' See \link{graph_prune} for details on the pruning procedure.
#'
#' @param m         A (square, possibly sparse) adjacency matrix.
#' @param vdf       An optional data frame with vertex metadata to add as vertex attributes.
#'                  Defaults to NULL.
#' @param mode      Character vector of length 1. Passed to
#'                  \code{igraph::graph_from_adjacency_matrix}. Defaults to 'undirected' (i.e. m
#'                  will be treated as symmetric).
#' @param weighted  Logical. Treat entries in m as edge weights. Defaults to TRUE.
#' @param normalize Logical. Normalize weights to the (0-1] range. Defaults to TRUE.
#' @param diag      Logical. Include diagonal entries. Passed to
#'                  \code{igraph::graph_from_adjacency_matrix}. Defaults to FALSE.
#'
#' @return An igraph object with the resulting graph.
#'
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix
graph_stitch <- function(
    m, vdf=NULL,
    mode='undirected', weighted=TRUE, normalize=TRUE, diag=FALSE, tol=1, verbose=FALSE
) {
    g <- m %>%
        igraph::graph_from_adjacency_matrix( mode=mode, weighted=weighted, diag=diag ) %>%
        graph_prune( tol=tol, verbose=verbose )
    if( normalize ) {
        if( !is.null( igraph::edge_attr( g, 'weight' ) ) ) {
            w <- igraph::edge_attr( g, 'weight' )
            igraph::edge_attr( g, 'weight' ) <- w / max( w )
        } else {
            warning( 'Can\'t normalize weights in unweighted graph!' )
        }
    }
    if( !is.null( vdf ) ) g %<>% graph_add_vertex_data( vdf )
    return( g )
}

#' Threshold pruning maintaining connectivity
#' 
#' This function will attempt to remove as many edges as possible following some order of edge 
#' significance, while maintaining total graph connectivity up to the given vertex dettachment 
#' tolerance.
#' 
#' Different prunning strategies can be implemented via different edge sorting orders according to
#' some edge significance criterion. Naive thresholding will consider each edge's weight as its 
#' significance.
#' 
#' Internally, this function will try to find an optimum significance threshold through a pivot 
#' search over the entire significance range. If \code{verbose} is \code{TRUE}, the progress of the 
#' pivot search is printed to output to show the algorithm's optimization process; this is 
#' sometimes useful to analyze degenerate results.
#'
#' NB: be advised that this strategy may fail for some graphs, e.g. if all vertices dettach from 
#' the main component in components smaller than the accepted tolerance, the algorithm will happily 
#' proceed to delete all edges in the graph. Lowering the connectivity threshold may sometimes help,
#' but in most cases a different edge sorting strategy is needed. See \link{graph_edge_score}.
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
graph_prune <- function(
    g, edges=NULL, tol=1, attr='weight', desc=TRUE, dropV=TRUE, verbose=FALSE
) {
    chk_igraph( g )

    if( is.null( edges ) ) {
      edges <- E( g )[ order( edge_attr( g, attr ), decreasing=desc ) ]
    }

    nv <- length( igraph::V( g ) )
    ne <- length( edges )
    if( verbose ) message( sprintf( 'Pruning graph to max dettached component size %d', tol ) )

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

    if( verbose ) message( sprintf(
        "Edge pruning removed %4.2f%% of edges and dropped %4.2f%% of vertices",
        ( ( ne - length( igraph::E( out ) ) ) / ne ) * 100,
        ( ( nv - length( igraph::V( out ) ) ) / nv ) * 100
    ) )

    return( out )
}

#' Determine graph connectivity up to the given tolerance
#'
#' This function will return true if the maximum component size of all components detached from the
#' largest component is greater or equal to the given tolerance. E.g. a tolerance of one implies
#' that a graph will be considered connected if all minor components are at most orphan nodes,
#' a tolerance of 2 implies the same if all minor components are at most only dyads, etc.
#'
#' This function will call \code{igraph::components( g )} if no value is given for \code{cmps}.
#' \emph{This can be extremely slow}, but so far I have not been able to find or produce
#' a faster component determination strategy.
#' If you know of one, please contact me (e.g. if you know how to compute matrix kernels in less
#' than \eqn{O(n^3)} )
#'
#' @param g    An igraph graph.
#' @param tol  An integer indicating the maximum size of detached components. Defaults to 1.
#' @param cmps An igraph components object. Will be computed over g if none given.
#'
#' @return \code{TRUE} if the given graph has no detached component larger than tol, \code{FALSE}
#'         otherwise.
#'
#' @export
graph_connected <- function( g, tol=1, cmps=igraph::components( g ) ) {
    if( cmps$no == 1 ) return( TRUE )
    if( cmps$no > 1 && tol == 0 ) return( false )
    if( max( cmps$csize[ -which.max( cmps$csize ) ] ) > tol ) {
        return( FALSE )
    } else {
        return( TRUE )
    }
}

#' Add vertex data from the given data frame
#'
#' This function will add all the values contained in the given data frame as vertex attributes,
#' retaining variable names as attribute names, optionally prefixed by the given \code{pref} string.
#'
#' Igraph does not accept all of \R's data types as attribute values, so some types are coerced
#' without information loss: logical values are converted to integers and factors are replaced with
#' their level names as strings. All other values are left unchanged.
#'
#' @param g    An igraph object.
#' @param vdf  A valid lexical data frame with vertex data. See \link{lexical_dataset}.
#' @param pref An optional prefix for attribute names.
#'
#' @return \code{g}, with data from vdf added as vertex attributes.
#'
#' @export
graph_add_vertex_data <- function( g, vdf, pref='' ) {
    data <- term_data( vdf, names( V( g ) ) )
    for( n in names( data ) ) {
        igraph::vertex_attr( g, pref %.% n ) <- if( is.logical( ( v <- data[[n]] ) ) ) {
            as.integer( v )
        } else if( is.factor( v ) ) {
            as.character( v )
        } else {
            v
        }
    }
    return( g )
}

#' Add cluster data to vertices and edges
#'
#' This function will add cluster data from the given \code{cms} igraph communities object as
#' vertex and edge attributes.
#'
#' Edge attributes added by this function consist of a logical value named 'xing' indicating wether
#' the respective edge crosses community boundaries or not, as well as a numeric "cluster weight"
#' value named the same as the given \code{eweight} name with a '_c' suffix, equal to the product
#' of the respective edge's value in the attribute \code{eweight} attribute and the value of
#' \code{intra_factor} if the edge is internal to a community and 1 otherwise.
#'
#' Vertex attributes added by this function consist of the vertex's community membership as an
#' integer value named 'comm'. If \code{contrib} is \code{TRUE}, then cluster-vertex contribution
#' scores will be computed and also added as numeric values named 'wgt_c2v' and 'wgt_v2c' for
#' cluster-vertex contribution and vertex-cluster contribution scores, respectively. NB: this
#' computation can be very slow for large networks.
#'
#' @param g            An igraph object.
#' @param cms          An igraph communities object.
#' @param intra_factor Factor to multiply edge weights for non-community crossing edges.
#'                     Defaults to 10.
#' @param eweight      Name for the input edge weight attribute. Defaults to 'weight'.
#' @param pref         An optional prefix that will be added to all attribute names.
#' @param contribs     Logical. Add cluster-vertex contribution scores. Defaults to \code{FALSE}
#' @param quiet        Logical. Suppress all messages. Defaults to \code{FALSE}.
#'
#' @return \code{g} with vertex and edge attributes containing relevant data from \code{cms}.
#'
#' @export
graph_add_cluster_data <- function(
    g, cms, intra_factor=10, eweight='weight', pref=NULL, contribs=FALSE, quiet=FALSE
) {
    pref <- if( is.null( pref ) || pref == '' ) '' else pref %.% '_'

    # edges
    x <- igraph::crossing( cms, g )
    igraph::edge_attr( g, pref %.% 'xing' ) <- x %>% as.integer()
    w <- igraph::edge_attr( g, eweight )
    igraph::edge_attr( g, pref %.% eweight %.% '_c' ) <- w * ifelse( x, intra_factor, 1 )

    # vertices
    k <- igraph::membership( cms )
    igraph::vertex_attr( g, pref %.% 'comm' ) <- k
    if( contribs ) {
        if( !quiet ) message( 'Computing cluster-vertex weights. This may a take a bit...' )
        igraph::vertex_attr( g, pref %.% 'wgt_v2c' ) <- graph_cluster_contribs( g, k, mode='vc' )
        igraph::vertex_attr( g, pref %.% 'wgt_c2v' ) <- graph_cluster_contribs( g, k, mode='cv' )
    }
    return( g )
}

#' Create signature for vertex sets.
#'
#' This function will create a unique value for each of the vertex sets defined in the given vector
#' k, by applying the given signature function to the vertex attributes contained in the given
#' graph to the set of topn elements in each set in k, filtered by fltr, sorted according to the
#' results of the given score function.
#'
#' By default, this function will filter by pos to select only nouns, will rank vertices according
#' to the product of their cluster contribution scores, and concantenate the term of the top 10
#' vertices in each cluster.
#'
#' @param g          An igraph graph
#' @param fltr       A vector to filter vertices in \code{g} before computing scores or sorting.
#'                   Defaults to \code{ pos == 'NN' | pos == 'NP'}, i.e. nouns.
#' @param k          A membership vector with sets to compute signatures for. Defaultf to 'comm',
#'                   i.e. communities.
#' @param score.func A function to compute vertex scores in each set in \code{k} from the
#'                   vertex attributes available in \code{g}. Defaults to the product of the
#'                   cluster-vertex contribution scores. See \link{graph_cluster_contribs}.
#' @param desc       Logical. Sort by the value of \code{score.func} in descending order. Defaults
#'                   to TRUE.
#' @param topn       Integer. The number of vertices to include in the signature. Defaults to 10.
#'                   Set to NA to include all vertices.
#' @param sig.func   A signature function that will receive the vertex attribute values for the
#'                   \code{topn} members in each set, and should return a character value to use as
#'                   set signature. Defaults to the concatenation of the vertices' terms.
#'
#' @export
graph_make_set_sigs <- function( g,
    fltr=( igraph::vertex_attr( g )$pos == 'NN' | igraph::vertex_attr( g )$pos == 'NP' ),
    k=igraph::vertex_attr( g, 'comm' ), score=function( df ) df$wgt_v2c * df$wgt_c2v, desc=TRUE,
    topn=10, sig.func=function( df ) df$term %>% paste( collapse=' ' )
) {
    df <- g %>% gr$vattr() %>% as.data.frame()
    df <- df[ fltr, ]
    df$crank <- vapply( 1:nrow( df ), function( r ) {
        score( df[ r, ] )
    }, 0.0 )
    topn <- ifelse( is.na( topn ), Inf, topn )
    sigs <- vapply( unique( k ), function( ki ) {
        rows <- df[ k[fltr] == ki, ][ order( df[ k[fltr] == ki, ]$crank, decreasing=desc ), ]
        rows[ 1:min( nrow( rows ), topn ), ] %>%
            sig.func() %>% return()
    }, '' )
    return( sigs )
}

# Clustering and cluster scores ---------------------------------------------------------------

#' Detect communities using the given clustering function.
#'
#' Constructs an igraph communities object by applying the given \code{clust_func} community
#' detection function over the given \code{g} igraph object.
#'
#' \code{clust_func} must follow igraph's clustering functions API, i.e. take a graph object as
#' input and produce an igraph 'communities' object as output.
#'
#' @param g            An igraph graph.
#' @param cluster_func A function implementing a clustering algorithm. Defaults to
#'                     \link{igraph::clustr_louvain}.
#' @param undirected   Logical. Ignore edge directionality, Defalts to TRUE.
#' @param quiet        Logical. Suppress progress messages. Defaults to FALSE.
#' @param ...          Additional parameters passed to \code{clust_func}.
#'
#' @return An igraph communities object with the clustering results.
#' 
#' @export
graph_cluster <- function(
    g, cluster.func=igraph::cluster_louvain, undirected=TRUE, quiet=FALSE, ...
) {
    if( undirected ) {
        g <- igraph::as.undirected( g )
    }
    s <- proc.time()
    c <- cluster.func( g, ... )
    e <- proc.time()
    if( !quiet ) message( sprintf(
        "Communities extracted in %8.4f seconds: %d groups, %6.4f mod.",
        ( e - s )[3], length( c ), igraph::modularity( c )
    ) )
    return( c )
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
#' @param g      An igraph graph or dense adjacency matrix.
#' @param k      A vertex-cluster membership vector.
#' @param mode   The direction of the measure; one of 'vc' or 'cv'. See details.
#' @param matrix Logical. Return vertex-cluster matrix instead of a score vector. Useful for
#'               cluster similarity computations. Defaults to \code{FALSE}
#'
#' @return If \code{matrix} is \code{FALSE}, a vector of length equal to the number of vertices in
#'         the given vset with the requested contribution measure values for each vertex in its
#'         assigned community.
#'         If \code{matrix} is \code{TRUE}, a matrix with as many rows as vertices in \code{vset}
#'         and as many columns as unique values in \code{k}.
#'
#' @export
graph_cluster_contribs <- function( ... ) {
    UseMethod( "graph_cluster_contribs" )
}

#' @export
graph_cluster_contribs.igraph <- function( g, ..., weight='weight' ) {
    m <- igraph::as_adj( g, type='both', attr=weight, sparse=FALSE )
    graph_cluster_contribs( m, ... )
}

#' @export
graph_cluster_contribs.matrix <- function( m, k, mode=c('cv','vc'), as.matrix=FALSE ) {
    mode  <- match.arg( mode )
    cwt <- switch( match.arg( mode ), cv=c2v_contrib( m, k ), vc=v2c_contrib( m, k ) )
    if( as.matrix ) {
        ks <- unique( k )
        out <- matrix( NA, nrow=length( V( g ) ), ncol=length( ks ) )
        for( i in 1:length( ks ) ) {
            out[,i] <- ifelse( k == ks[i], cwt, 0 )
        }
    } else {
        out <- cwt
    }
    return( out )
}

# Edge scoring --------------------------------------------------------------------------------

#' Compute edge significance scores with the given function
#' 
#' This function will apply \code{score.func} over all edges in a graph, passing as parameters each 
#' edge's weight, the total incoming strength of its head vertex, the total outgoing strength of 
#' its tail vertex and the total weight in the graph, for e.g. compuation of some probabilistic 
#' edge significance measure. See \link{graph_edge_score_ident} for a null implementation that 
#' simply returns the weight; See \link{graph_edge_score_mlf} for an implemenation of Dianati's
#' Marginal Likelihood Filter.
#' 
#' TODO: fast matrix -> igraph edge attribute mappings would allow for considerable speed-up.
#' 
#' @param g          An igraph graph.
#' @param score.func An edge scoring function with signature
#'                   \code{function( w, head_s, tail_S, total_w )} giving edge's weight, total head 
#'                   vertex strength, total tail vertex strength and total graph weigth.
#' @param src.attr   An attribute to read edge weights from. Defaults to \code{'weight'}.
#' @param tgt.attr   A target attribute to write scores to. Defaults to \code{NULL}: don't write 
#'                   scores back to graph.
#' 
#' @return A vector of edge scores.
#' 
#' @export
graph_edge_score <- function(
    g, score.func=graph_edge_score_mlf, src.attr='weight', tgt.attr=NULL
) {
    edges   <- E( g )
    heads   <- head_of( g, edges )
    tails   <- tail_of( g, edges )
    head_s  <- strength( g, heads, mode='in' )
    tail_s  <- strength( g, heads, mode='out' )
    weights <- edge_attr( g, src.attr )
    browser()
    res <- score.func( weights, head_s, tail_s, sum( weights ) )
    if( !is.null( tgt.attr ) ) {
        edge_attr( g, tgt.attr ) <- res
    }
    return( res )
}

#' No-op edge scoring function.
#' 
#' Does nothing: score edges according to their given weight \code{w}.
#' 
#' @param w       A vector of edge weights. Returned as is.
#' @param head_s  A vector of edge's head vertices' total strength. Ignored.
#' @param tail_s  A vector of edge's tail vertices' total strength. Ignored.
#' @param total_w The total weight of edges in the source graph. Ignored.
#' 
#' @return This is a no-op function that returns \code{w}.
#' 
#' @export
graph_edge_score_ident <- function( w, head_s, tail_s, total ) {
    return( w )
}

#' Navid Dianati's Marginal Likelihood Filter.
#' 
#' Computes edge significances as the p-value associated to the given edges weight as the 
#' outcome of a set of \code{total_w} Bernoulli trials with probability equal to 
#' \eqn{head_s * tail_s / total_w ^2}.
#' 
#' Since Dianati's original implementation was designed for integer-weighted edges (i.e. 
#' multigraphs modeling actual ) using the binominal distribution, this function will naively 
#' convert all parameters to integers before computing the p-value.
#' 
#' TODO: Citation needed.
#' 
#' @param w       A vector of edge weights.
#' @param head_s  A vector of edge's head vertices' total strength.
#' @param tail_s  A vector of edge's tail vertices' total strength.
#' @param total_w The total weight of edges in the source graph.
#' 
#' @export
graph_edge_score_mlf <- function( w, head_s, tail_s, total_w ) {
    x <- w %>% floor() %>% as.integer()
    n <- total %>% floor() %>% as.integer()
    p <- ( head_s * tail_s ) / ( total^2 )
    r <- binom.test( x, n, p, alternative="greater" )
    return( r )
}


# Cluster alignment ---------------------------------------------------------------------------

#' Cluster alignment
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
        if( is.null( igraph::vertex_attr( g1, vid ) ) || is.null( igraph::vertex_attr( g2, vid ) ) ) {
            stop( "Invalid name for vertex id attr, needed for cross-graph alignment" )
        }
        if( fill ) { # union
            unique( c( igraph::vertex_attr( g1, vid ), igraph::vertex_attr( g2, vid ) ) )
        } else { # intersection
            igraph::vertex_attr( g1, vid ) %in% igraph::vertex_attr( g2, vid )
        }
    } else {
        igraph::V( g1  ) %>% as.integer()
    }

    g2 <- if( is.null( g2 ) ) g1 else g2
    m1 <- contrib_func( g1, cms1, ... )
    m2 <- contrib_func( g2, cms2, ... )
    rownames( m1 ) <- igraph::vertex_attr( g1, vid )
    rownames( m2 ) <- igraph::vertex_attr( g2, vid )

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
    gm[ 1:k1, ( k1 + 1 ):ncol( gm ) ] <- m
    gm[ ( k1 + 1 ):ncol( gm ), 1:k1 ] <- t( m )
    g <- igraph::graph_from_adjacency_matrix( gm, weighted=TRUE )
    igraph::vertex_attr( g, "type" ) <- c( rep( TRUE, k1 ), rep( FALSE, k2 ) )
    return( g )
}

# Non-export utility functions ----------------------------------------------------------------

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
