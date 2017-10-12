#include "innerp.hpp"

#if INNERP_IMPL == IMPL_SERIAL
#include "tools.hpp"

double innerp( Vec vi, Vec vj, ind i, ind j, int mode ) {
    return get_func( mode )( vi, vj, i, j );
}

Vec marg_innerp( Mat a, Mat b, int mode ) {
    if( a.rows() != b.rows() || a.cols() != b.cols() ) {
        Rcpp::stop( "wrong dimensions in input: a and b must have the same shape" );
    }
    Vec out = Vec::Zero( a.rows() );
    for( ind i = 0; i < a.rows(); i++ ) {
        out[ i ] = innerp( a.row( i ), b.row( i ), i, i, mode );
    }
    return out;
}

Mat full_innerp( Mat a, Mat b, int mode ) {
    if( a.cols() != b.cols() ) {
        Rcpp::stop( "wrong dimensions in input: a cols and b rows must be equal" );
    }
    Mat out = Mat::Zero( a.rows(), b.rows() );
    for( ind i = 0; i < a.rows(); i++ ) {
        for( ind j = 0; j < b.rows(); j++ ) {
            out( i, j ) = innerp( a.row( i ), b.row( j ), i, j, mode );
        }
    }
    return out;
}

Mat self_innerp( Mat a, int mode, bool self ) {
    Mat out = Mat::Zero( a.rows(), a.rows() );
    bool symm = symmetric( mode );
    for( ind i = 0; i < a.rows(); i++ ) {
        for( ind j = 0; j < ( symm ? i + ( self ? 1 : 0 ) : a.rows() ); j++ ) {
            if( !self && i == j ) continue; // needed for !self && !symm
            double v = get_func( mode )( a.row( i ), a.row( j ), i, j );
            out( i, j ) = v;
            if( symm ) out( j, i ) = v;
        }
    }
    return out;
}

#endif
