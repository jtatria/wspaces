#ifndef SIMDIV_H
#define SIMDIV_H 1

#include "wspaces_types.hpp"
#include "tools.hpp"

typedef F<Vec,Vec,ind,ind > SimFunc;

enum Sim {
    ADDITIVE  = 0,
    DWEIGHTED = 1,
    COSINE    = 2,
    KL_DIV    = 3,
    JS_DIV    = 4,
    BHATTA    = 5,
    HELLINGER = 6,
    JACCARD   = 7,
};

//' @export
// [[Rcpp::export]]
double sim_additive( const Vec& w1, const Vec& w2, ind i, ind j ) {
    double num = ( w1.array() > 0 && w2.array() > 0 ).select( ( w1 + w2 ), 0 ).sum();
    double den = w1.sum();
    return num / den;
}

//' @export
// [[Rcpp::export]]
double sim_dweighted( const Vec& w1, const Vec& w2, ind i, ind j ) {
    double num = ( w1.cwiseMin( w2 ).array().sum() );
    double den = w1.sum();
    return num / den;
}

//' Cosine similarity.
//'
//' The cosine similarity is equal to the vector dot product between v and u divided by the product
//' of their norms: \eqn{ cos( v, u ) = \frac{}v \cdot u}{\lVert v \rVert \lVert u \rVert } }.
//'
//' Note that this measure is symmetric and induces a proper distance when taking its recyprocal.
//'
//' @param v A vector.
//' @param u A vector.
//'
//' @return A scalar value \eqn{ cos( v, u ) \in [ -1, 1 ] }.
//'
//' @export
// [[Rcpp::export]]
double sim_cosine( const Vec& v, const Vec& u ) {
    if( v.isZero() || u.isZero() ) return 0;
    double num = v.dot( u );
    double den = v.norm() * u.norm();
    return num / den;
}

//' Kullback-Leibler divergence
//'
//' The Kullback-Leibler divergence of p from q is equal to
//' \eqn{ \sum_{i} p_i\ log\ { \frac{p_i}{q_i} } }.
//' KL is valid only when q_i = 0 iff p_i = 0. This function will fail if p or q contain any zeros.
//'
//' @param p     A frequency or probability vector. Coerced to prob by p / sum( p )
//' @param q     A frequency or probability vector. Coerced to prob by q / sum( q )
//'
//' @return A scalar value equal to the Kullback-Leibler divergence of p from q.
//'
//' @export
// [[Rcpp::export]]
double div_kulback_leibler( const Vec& p, const Vec& q, bool coerce=false ) {
    Vec pi = p.array() / p.sum();
    Vec pj = q.array() / q.sum();
    return ( ( pi.array() / pj.array() ).log() * pi.array() ).sum();
}

//' Jensen-Shannon divergence
//'
//' Computes the Jensen-Shannon smoothing reflection of the Kullback-Leibler divergence between the
//' probability distributions represented by the given frequency or probabality vectors. The
//' Jensen-Shannon divergence between p and q is equal to the average of the Kullback-Leibler
//' divergence between p and m and between q and m, where m is the point-wise average between p and
//' q. i.e.: \eqn{ JS_{pq} = ( KL_{p, (p + q) / 2} + KL_{q, (p + q) / 2}) / 2 }.
//'
//' This function only computes m and the average between the two Kullback-Leibler divergences
//' between p and q and m. The Kullback-Leibler divergences are computed by
//' \code{div_kullback_leibler}.
//'
//' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
//' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
//'
//' @return a scalar value equal to the Jensen-Shannon divergence between p and q.
//'
//' @export
// [[Rcpp::export]]
double div_jensen_shannon( const Vec& p, const Vec& q, bool coerce=false ) {
    Vec m = ( p + q ) / 2.0;
    double kl_p = div_kulback_leibler( p, m );
    double kl_q = div_kulback_leibler( q, m );
    return( ( kl_p + kl_q ) / 2.0 );
}

//' Bhattacharyya coefficient
//'
//' Computes the Bhattacharyya overlap coefficient between the probability distributions represented
//' by the given frequency or probability vectors. The Bhattacharyya coefficient is equal to
//' \eqn{\sum_{i}^{n} \sqrt{ p_i*q_i } }.
//'
//' The given vectors will be coerced to probabilities by \code{x = x / sum( x )}.
//'
//' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
//' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
//'
//' @return a scalar value equal to the Bhattacharyya overlap coefficient between p and q.
//'
//' @export
// [[Rcpp::export]]
double coef_bhattacharyya( const Vec& p, const Vec& q ) {
    Vec pi = p.array() / p.sum();
    Vec pj = q.array() / q.sum();
    return ( ( pi.array() * pj.array() ).sqrt() ).sum();
}

//' Bhattacharyya divergence
//'
//' Computes the Bhattacharyya divergence between the probability distributions represented by the
//' given frequency or probabality vectors. The Bhattacharyya divergence is equal to the reciprocal
//' of the natural logarithm of the Bhattacharyya coefficient between p and q:
//' \eqn{ -log( \sum_{i}^{n} \sqrt{ p_i*q_i } ) }
//'
//' This function just computes the natural logarithm reciprocal. The coefficient is computed by
//' \code{coef_bhattacharyya( p, q )}.
//'
//' This measure is not a proper distance because it does not respect the triangle inequality.
//' If a proper distance is needed, use the Hellinger distance (\code{dist_hellinger( p, q )}), also
//' based on the Bhattacharyya coefficient.
//'
//' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
//' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
//'
//' @return a scalar value equal to the Bhattacharyya divergence between p and q.
//'
//' @export
// [[Rcpp::export]]
double div_bhattacharyya( const Vec& p, const Vec& q ) {
    double coef = coef_bhattacharyya( p, q );
    return std::log<double>( coef ).real() * -1;
}

//' Hellinger distance
//'
//' Computes the Hellinger distance between the probability distributions represented by the given
//' frequency or probability vectors p and q. The Hellinger distance is equal to the square root of
//' the 1-complement of the Bhattacharyya coefficient between p and q:
//' \eqn{ \sqrt{ 1 - \sum_{i}^{n} \sqrt{ p_i*q_i } } }.
//'
//' This function only computes the square root of the 1-complement. The coefficient is computed by
//' \code{coef_bhattacharyya(p, q)}.
//'
//' Unlike the closely related Bhattacharyya divergence, this measure is a proper distance as it
//' respects the triangle inequality.
//'
//' @param p A frequency or probability vector. Coerced to prob by p / sum( p )
//' @param q A frequency or probability vector. Coerced to prob by q / sum( q )
//'
//' @return a scalar value equal to the Hellinger distance between p and q.
//'
//' @export
// [[Rcpp::export]]
double dist_hellinger( const Vec& p, const Vec& q ) {
    double bha = coef_bhattacharyya( p, q );
    return std::sqrt<double>( 1 - bha ).real();
}

//' @export
// [[Rcpp::export]]
double dist_jaccard( const Vec& a, const Vec& b ) {
    double num = ( a.cwiseMin( b ) ).sum();
    double den = ( a.cwiseMax( b ) ).sum();
    return( num / den );
}

inline bool symmetric( const Sim sim ) {
    switch( sim ) {
        case ADDITIVE  : return false;
        case DWEIGHTED : return false;
        case COSINE    : return true;
        case KL_DIV    : return false;
        case JS_DIV    : return true;
        case BHATTA    : return true;
        case HELLINGER : return true;
        case JACCARD   : return true;
        default : Rcpp::stop( "Unknown sim/div function requested" );
    }
}

//' @export
// [Rcpp::export]]
bool simdiv_symm( int mode ) {
    return( symmetric( static_cast<Sim>( mode ) ) );
}

inline SimFunc get_func( const Sim sim ) {
    switch( sim ) {
        case ADDITIVE  : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return sim_additive( vi, vj, i, j );
        };
        case DWEIGHTED : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return sim_dweighted( vi, vj, i, j );
        };
        case COSINE    : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return sim_cosine( vi, vj );
        };
        case KL_DIV    : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return div_kulback_leibler( vi, vj );
        };
        case JS_DIV    : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return div_jensen_shannon( vi, vj );
        };
        case BHATTA    : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return div_bhattacharyya( vi, vj );
        };
        case HELLINGER : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return dist_hellinger( vi, vj );
        };
        case JACCARD   : return []( const Vec& vi, const Vec& vj, ind i, ind j ) -> double {
            return dist_jaccard( vi, vj );
        };
        default : Rcpp::stop( "Unknown sim/div function requested" );
    }
}

#endif
