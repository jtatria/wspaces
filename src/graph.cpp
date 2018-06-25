#include "tools.hpp"
#include <RcppParallel.h>
#include <unsupported/Eigen/KroneckerProduct>

namespace impl {
    struct ContribWorker : public RcppParallel::Worker {
        RcppParallel::RVector<double> tgt;
        F<ind> func;

        ContribWorker( RVecD tgt, F<ind> func ) : tgt( tgt ), func( func ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( ind i = begin; i < end; ++i ) {
                tgt[i] = func( i );
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
