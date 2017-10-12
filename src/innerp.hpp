#ifndef INNERP_H
#define INNERP_H 1

#include "wspaces_types.hpp"

#define IMPL_SERIAL   0
#define IMPL_PARALLEL 1
#define IMPL_CUDA     2

// Change implementation compilation here
#define INNERP_IMPL IMPL_PARALLEL

typedef F<Vec,Vec,ind,ind > InnerP_func;

enum InnerP {
    DOT    = 0,
    COSINE = 1,
};

InnerP_func func_dot();

InnerP_func func_cosine();

inline InnerP_func get_func( int mode ) {
    switch( static_cast<InnerP>( mode ) ) {
        case DOT    : return func_dot();
        case COSINE : return func_cosine();
        default : Rcpp::stop( "Unknown inner product requested" );
    }
}

inline bool symmetric( int mode ) {
    switch( static_cast<InnerP>( mode ) ) {
        case DOT    : return true;
        case COSINE : return true;
        default : Rcpp::stop( "Unknown inner product requested" );
    }
}

// [[Rcpp::export]]
double innerp( Vec vi, Vec vj, ind i, ind j, int mode = 0 );

// [[Rcpp::export]]
Vec marg_innerp( Mat a, Mat b, int mode = 0 );

// [[Rcpp::export]]
Mat full_innerp( Mat a, Mat b, int mode = 0 );

// [[Rcpp::export]]
Mat self_innerp( Mat a, int mode = 0, bool self = true );

#endif



