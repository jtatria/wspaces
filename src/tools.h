#ifndef TOOLS_
#define TOOLS_ 1

// [[Rcpp::plugins("cpp11")]]
#include "wspaces_types.h" // Needed for RcppExports compilation.

const double EPSILON = 0.00000001;

template <typename... args >
using F = std::function<double(args...) >;

enum Margin {
    Col = 0,
    Row = 1,
};

Vec marg_sum( const SpMat&, const Margin& );
Vec marg_avg( const SpMat&, const Margin& );
Vec marg_prb( const SpMat&, const Margin& );
Vec marginal( const SpMat&, const Margin&, const F<double,double> );

#endif // TOOLS_
