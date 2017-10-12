#include "innerp.hpp"

#if INNERP_IMPL == IMPL_PARALLEL
#include "tools.hpp"
#include <RcppParallel.h>
#include <iostream>

struct MargWrkr : public RcppParallel::Worker {
    const Mat& a;
    const Mat& b;
    const InnerP_func func;
    Vec& tgt;

    MargWrkr( Mat& a, Mat& b, Vec& tgt, InnerP_func func ) :
        a( a ), b( b ), tgt( tgt ), func( func ) {
    }

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            tgt[i] = func( a.row( i ), b.row( i ), i, i );
        }
    }
};

struct FullWrkr : public RcppParallel::Worker {
    const Mat& a;
    const Mat& b;
    const InnerP_func func;
    Mat& tgt;

    FullWrkr( Mat& a, Mat& b, Mat& tgt, InnerP_func func ) :
        a( a ), b( b ), tgt( tgt ), func( func ) {
        Rcpp::Rcout << "a rows: " << a.rows() << std::endl;
        Rcpp::Rcout << "b rows: " << b.rows() << std::endl;
        Rcpp::Rcout << "tgt dim: " << tgt.rows() << " " << tgt.cols() << std::endl;
    }

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            for( ind j = 0; j < b.rows(); j++ ) {
                double v = func( a.row( i ), b.row( j ), i, j );
                tgt.coeffRef( i, j ) = v;
            }
        }
    }
};

struct SelfWrkr : public RcppParallel::Worker {
    const Mat& a;
    const InnerP_func func;
    const bool symm;
    const bool self;
    Mat& tgt;

    SelfWrkr( Mat& a, Mat& tgt, InnerP_func func, bool symm, bool self ) :
        a( a ), tgt( tgt ), func( func ), symm( symm ), self( self ) {
    }

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            for( ind j = 0; j < ( symm ? i + ( self ? 1 : 0 ) : a.rows() ); j++ ) {
                if( !self && i == j ) continue; // needed for !self && !symm
                double v = func( a.row( i ), a.row( j ), i, j );
                tgt.coeffRef( i, j ) = v;
                if( symm ) tgt.coeffRef( j, i ) = v;
            }
        }
    }
};

double innerp( Vec vi, Vec vj, ind i, ind j, int mode ) {
    return get_func( mode )( vi, vj, i, j );
}

Vec marg_innerp( Mat a, Mat b, int mode ) {
    if( a.rows() != b.rows() || a.cols() != b.cols() ) {
        Rcpp::stop( "wrong dimensions in input: a and b must have the same shape" );
    }
    Vec out = Vec::Zero( a.rows() );
    MargWrkr wrkr( a, b, out, get_func( mode ) );
    RcppParallel::parallelFor( 0, a.rows(), wrkr );
    return out;
}

Mat full_innerp( Mat a, Mat b, int mode ) {
    if( a.cols() != b.cols() ) {
        Rcpp::stop( "wrong dimensions in input: a cols and b rows must be equal" );
    }
    Mat out = Mat::Zero( a.rows(), b.rows() );
    FullWrkr wrkr( a, b, out, get_func( mode ) );
    RcppParallel::parallelFor( 0, a.rows(), wrkr );
    return out;
}

Mat self_innerp( Mat a, int mode, bool self ) {
    Mat out = Mat::Zero( a.rows(), a.rows() );
    SelfWrkr wrkr( a, out, get_func( mode ), symmetric( mode ), self );
    RcppParallel::parallelFor( 0, a.rows(), wrkr );
    return out;
}

#endif
