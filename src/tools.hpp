#ifndef TOOLS_
#define TOOLS_ 1

#include "wspaces_types.hpp" // Needed for RcppExports compilation.

const double EPSILON = 0.00000001;

template <typename... args >
using F = std::function<double(args...)>;

enum Margin {
    Col = 0,
    Row = 1,
};

Vec marg_sum( const SpMat&, const Margin& );
Vec marg_avg( const SpMat&, const Margin& );
Vec marg_prb( const SpMat&, const Margin& );
Vec marginal( const SpMat&, const Margin&, const F<double,double> );

inline bool isomorphic( const Mat m1, const Mat m2, bool bail=true ) {
    if( m1.rows() == m2.rows() && m1.cols() == m2.cols() ) {
        return true;
    } else {
        if( bail ) Rcpp::stop( "Objects are not isomorphic" );
        else return false;
    }
}

#endif // TOOLS_
