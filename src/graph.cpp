#include "tools.hpp"
#include <RcppParallel.h>
#include <unsupported/Eigen/KroneckerProduct>

namespace impl {
    struct ContribWorker : public RcppParallel::Worker {
        RcppParallel::RVector<double> m_tgt;
        F<ind> m_func;

        ContribWorker( RVecD tgt, F<ind> func ) : m_tgt( tgt ), m_func( func ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( ind i = begin; i < end; ++i ) {
                m_tgt[i] = m_func( i );
            }
        }
    };

    struct EdgeScoreWorker : public RcppParallel::Worker {
        RcppParallel::RMatrix<double> m_tgt;
        F<ind,ind> m_func;

        EdgeScoreWorker( RMatD tgt, F<ind,ind> func ) : m_tgt( tgt ), m_func( func ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( ind i = begin; i < end; ++i ) {
                for( ind j = 0; j < m_tgt.nrow(); ++j ) {
                    m_tgt.row( i )[j] = m_func( i, j );
                }
            }
        }
    };
};

//' @export
// [[Rcpp::export]]
RVecD c2v_contrib( const Mat& adj, const Ivec& k_memb ) {
    if( adj.rows() != adj.cols() ) Rcpp::stop( "adj is not a square!" );
    Ivec no   = Ivec::Zero( k_memb.size() );
    Ivec yes  = Ivec::Ones( k_memb.size() );
    const Imat zero = Imat::Zero( adj.rows(), adj.cols() );
    F<ind> func = [&]( ind i ) -> scalar {
        int k = k_memb( i );
        scalar num_in  = ( k_memb.array() == k ).select( adj.row( i ), no ).sum();
        scalar num_out = ( k_memb.array() == k ).select( adj.col( i ), no ).sum();
        // all internal edges in k
        Ivec k_vs = ( k_memb.array() == k ).select( yes, no );
        Imat mask = Eigen::kroneckerProduct( k_vs, k_vs.transpose() ).cwiseSign();
        scalar den = ( mask.array() == 1 ).select( adj, zero ).sum();
        return ( num_in + num_out ) / ( den * 2 );
    };
    RVecD r( adj.rows() );
    impl::ContribWorker wrkr( r, func );
    RcppParallel::parallelFor( 0, adj.rows(), wrkr );
    return( r );
}


//' @export
// [[Rcpp::export]]
RVecD v2c_contrib( const Mat& adj, const Ivec& k_memb ) {
    if( adj.rows() != adj.cols() ) Rcpp::stop( "adj is not a square!" );
    Ivec no   = Ivec::Zero( k_memb.size() );
    F<ind> func = [&]( ind i ) -> scalar {
        int k = k_memb( i );
        // incident internal to k
        scalar num_in  = ( k_memb.array() == k ).select( adj.row( i ), no ).sum();
        scalar num_out = ( k_memb.array() == k ).select( adj.col( i ), no ).sum();
        // all incident
        scalar den = adj.row( i ).sum() + adj.col( i ).sum();
        return ( num_in + num_out ) / den;
    };
    RVecD r( adj.rows() );
    impl::ContribWorker wrkr( r, func );
    RcppParallel::parallelFor( 0, adj.rows(), wrkr );
    return( r );
}

RMatD edge_score_mlf( const Mat& adj, bool directed=false ) {
    F<ind,ind> func;
    if( directed ) {
        func = [&]( ind i, ind j ) -> scalar {
            return 0.0; // TODO
        };
    } else {
        func = [&]( ind i, ind j ) -> scalar {
            return 0.0; // TODO
        };
    }
    RMatD r( adj.rows(), adj.cols() );
    impl::EdgeScoreWorker wrkr( r, func );
    RcppParallel::parallelFor( 0, adj.rows(), wrkr );
    return r;
}
