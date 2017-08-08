#include <RcppEigen.h>
#include <RcppParallel.h>

using namespace Rcpp;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

// [[Rcpp::export]]
S4 serial_test( S4 m ) {
    SpMat src = as<MSpMat>( m );
    SpMat tgt = as<MSpMat>( clone( m ) );
    for( int i = 0; i < src.rows(); ++i ) {
        for( SpMat::InnerIterator srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
            tgtIt.valueRef() = -1;
        }
    }
    S4 out = wrap( tgt );
    return out;
}

struct MyWorker : public RcppParallel::Worker {
    SpMat src;
    SpMat tgt;

    MyWorker( const SpMat& src, SpMat& tgt )
        : src( src ), tgt( tgt ) {}

    void operator()( std::size_t begin, std::size_t end ) {
        for( int i = begin; i < end; ++i ) {
            for( SpMat::InnerIterator srcIt( src, i ), tgtIt( tgt, i ); srcIt; ++srcIt, ++tgtIt ) {
                tgtIt.valueRef() = -1;
            }
        }
    }
};

// [[Rcpp::export]]
S4 para_test( S4 m ) {
    SpMat src = as<MSpMat>( m );
    SpMat tgt = as<MSpMat>( clone( m ) );
    MyWorker wrkr( src, tgt );
    RcppParallel::parallelFor( 0, src.outerSize(), wrkr );
    S4 out = wrap( tgt );
    return out;
}
