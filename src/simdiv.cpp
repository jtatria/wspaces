/*
 * Copyright (C) 2017 José Tomás Atria <jtatria at gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/***************************************************************************************************
      Implementation of similarity and divergence measures on dense matrices and vectors
***************************************************************************************************/

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends( RcppParallel )]]
#include "tools.hpp"
#include <functional>
#include <RcppParallel.h>

typedef F<Vec,Vec,ind,ind > SimFunc;

using namespace Rcpp;

enum Sim {
    ADDITIVE  = 0,
    DWEIGHTED = 1,
    COSINE    = 2,
    KL_DIV    = 3,
    JS_DIV    = 4,
    BHATTA    = 5,
    HELLINGER = 6,
};

inline Mat check_input( Mat m, Sim sim ) {
    switch( sim ) {
        case ADDITIVE: case DWEIGHTED: case COSINE: case BHATTA: case HELLINGER: {
            return m;
        }
        case KL_DIV: case JS_DIV: {
            if( ( m.array() == 0 ).any() ) {
                // Dirichlet prior.
                return ( m.array() + 1 ) / ( m.sum() + ( m.rows() * m.cols() ) );
            } else {
                return m;
            }
        }
    }
}

// [[Rcpp::export]]
double sim_dweighted( Vec vi, Vec vj, ind i, ind j ) {
    isomorphic( vi, vj );
    double ii = vi( i ), ij = vi( j ), ji = vj( i ), jj = vj( j );
    double mi = std::min<double>( ii, ji );
    double mj = std::min<double>( ji, jj );
    double num = ( vi.cwiseMin( vj ).array().sum() ) - ( mi + mj );
    double den = vi.sum() - ( ii + ij );
    return num / den;
}

// [[Rcpp::export]]
double sim_additive( Vec vi, Vec vj, ind i, ind j ) {
    isomorphic( vi, vj );
    double ii = vi( i ), ij = vi( j ), ji = vj( i ), jj = vj( j );
    double num = ( vi.array() > 0 && vj.array() > 0 ).select( vi, 0 ).sum() - ( ii + ij + ji + jj );
    double den = vi.sum() - ( ii + ij );
    return num / den;
}

// [[Rcpp::export]]
double sim_cosine( Vec vi, Vec vj ) {
    isomorphic( vi, vj );
    double num = vi.dot( vj );
    double den = vi.norm() * vj.norm();
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
// [[Rcpp::export]]
double div_kulback_leibler( Vec vi, Vec vj, bool coerce=false ) {
    isomorphic( vi, vj );
    Vec pi = vi.array() / vi.sum();
    Vec pj = vj.array() / vj.sum();
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
// [[Rcpp::export]]
double div_jensen_shannon( Vec vi, Vec vj, bool coerce=false ) {
    isomorphic( vi, vj );
    Vec m = ( vi + vj ) / 2.0;
    double kl_vi = div_kulback_leibler( vi, m );
    double kl_vj = div_kulback_leibler( vj, m );
    return( ( kl_vi + kl_vj ) / 2.0 );
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
// [[Rcpp::export]]
double coef_bhattacharyya( Vec vi, Vec vj ) {
    isomorphic( vi, vj );
    Vec pi = vi.array() / vi.sum();
    Vec pj = vj.array() / vj.sum();
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
// [[Rcpp::export]]
double div_bhattacharyya( Vec vi, Vec vj ) {
    isomorphic( vi, vj );
    double coef = coef_bhattacharyya( vi, vj );
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
//' @return a scalar value equal to the Hellinger divergence between p and q.
//'
// [[Rcpp::export]]
double dist_hellinger( Vec vi, Vec vj ) {
    isomorphic( vi, vj );
    double bha = coef_bhattacharyya( vi, vj );
    return std::sqrt<double>( 1 - bha ).real();
}

struct IntraWrkr : public RcppParallel::Worker {
    const Mat src;
    const bool symm;
    const bool self;
    const SimFunc func;

    RcppParallel::RMatrix<double> tgt;

    IntraWrkr( const Mat src, RMatD tgt, bool self, bool symm, SimFunc func )
        : src( src ), tgt( tgt ), func( func ), self( self ), symm( symm ) {
    }

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            for( std::size_t j = 0; j < ( symm ? i + ( self ? 1 : 0 ) : src.rows() ); j++ ) {
                if( !self && i == j ) continue; // needed for !self && !symm
                double v = func( src.row( i ), src.row( j ), i, j );
                tgt( i, j ) = v;
                if( symm ) tgt( j, i ) = v;
            }
        }
    }
};

struct InterWrkr : public RcppParallel::Worker {
    const Mat mi;
    const Mat mj;
    const SimFunc func;

    RcppParallel::RVector<double> tgt;

    InterWrkr( Mat mi, Mat mj, RVecD tgt, SimFunc func )
        : mi( mi ), mj( mj ), tgt( tgt ), func( func ) {
    }

    void operator()( std::size_t begin, std::size_t end ) {
        for( std::size_t i = begin; i < end; i++ ) {
            tgt[i] = func( mi.row( i ), mj.row( i ), i, i );
        }
    }
};

inline SimFunc get_func( const Sim sim ) {
    switch( sim ) {
        case ADDITIVE  : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return sim_additive( vi, vj, i, j );
        };
        case DWEIGHTED : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return sim_dweighted( vi, vj, i, j );
        };
        case COSINE    : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return sim_cosine( vi, vj );
        };
        case KL_DIV    : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return div_kulback_leibler( vi, vj );
        };
        case JS_DIV    : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return div_jensen_shannon( vi, vj );
        };
        case BHATTA    : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return div_bhattacharyya( vi, vj );
        };
        case HELLINGER : return []( Vec vi, Vec vj, ind i, ind j ) -> double {
            return dist_hellinger( vi, vj );
        };
        default : stop( "Unknown sim/div function requested" );
    }
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
        default : stop( "Unknown sim/div function requested" );
    }
}

//' Compute self similarity or divergence/distance measures on dense matrices.
//'
//' This function will compute the designated similarity or distance function between all row or
//' column vectors in the given matrix. The implementation assumes observations are stored as rows,
//' with features stored in columns, as per the usual R data frame format; setting transpose to
//' TRUE will interpret the given matrix in the opposite direction.
//'
//' Available measures:
//' \itemize{
//'     \item{0: Additive cooccurrence retrieval. See \code{sim_additive}.}
//'     \item{1: Difference-weighted cooccurrence retrieval. See \code{sim_dweighted}.}
//'     \item{2: Cosine similarity. See \code{sim_cosine}.}
//'     \item{3: Kullback-Leibler divergence. See \code{div_kullback_leibler}.}
//'     \item{4: Jensen-Shannon divergence. See \code{div_jensen_shannon}.}
//'     \item{5: Bhattacharyya divergence. See \code{div_bhattacharyya}.}
//'     \item{6: Hellinger distance. See \code{dist_hellinger}.}
//' }
//' @param m         A matrix representing a vector set to compute measures over.
//' @param self      Logical indicating whether self measures should be computed too.
//' @param mode      The desired measure. See details.
//' @param transpose Logical indicating whether the given matrix should be transposed before
//'                  computation. Defaults to FALSE (i.e. compute measure over row vectors).
//'
//' @return A square matrix of dimension equal to the number of observation vectors
//'         (rows by default) containing the values of the desired measure between all observations
//'         in the given matrix.
//'
// [[Rcpp::export]]
RMatD self_simdiv( Mat m, bool self=false, int mode=2, bool transpose=false ) {
    Sim sim = static_cast<Sim>( mode );
    m = transpose ? m.transpose() : m;
    m = check_input( m, sim );
    bool symm = symmetric( sim );
    SimFunc func = get_func( sim );
    RMatD out( m.rows(), m.rows() );
    IntraWrkr wrkr( m, out, self, symm, func );
    RcppParallel::parallelFor( 0, m.rows(), wrkr );
    return out;
}

//' Compute row or column wise similarity or divergence/distance measures across two samples.
//'
//' This function will compute the designated similarity or distance function between the
//' corresponding row or column vectors in the given sample matrices.
//' The implementation assumes observations are stored as rows,
//' with features stored in columns, as per the usual R data frame format; setting transpose to
//' TRUE will interpret the given matrix in the opposite direction.
//'
//' Available measures:
//' \itemize{
//'     \item{0: Additive cooccurrence retrieval. See \code{sim_additive}.}
//'     \item{1: Difference-weighted cooccurrence retrieval. See \code{sim_dweighted}.}
//'     \item{2: Cosine similarity. See \code{sim_cosine}.}
//'     \item{3: Kullback-Leibler divergence. See \code{div_kullback_leibler}.}
//'     \item{4: Jensen-Shannon divergence. See \code{div_jensen_shannon}.}
//'     \item{5: Bhattacharyya divergence. See \code{div_bhattacharyya}.}
//'     \item{6: Hellinger distance. See \code{dist_hellinger}.}
//' }
//' @param mi        A matrix representing a vector set for sample i.
//' @param mj        A matrix representing a vector set for sample i.
//' @param mode      The desired measure. See details.
//' @param transpose Logical indicating whether the given matrices should be transposed before
//'                  computation. Defaults to FALSE (i.e. compute measure over row vectors).
//'
//' @return A vector of dimension equal to the number of observations in the given matrices
//'         (rows by default) containing the values of the desired measure for all observations
//'         in the two samples.
//'
// [[Rcpp::export]]
RVecD cross_simdiv( Mat mi, Mat mj, int mode=2 ) {
    isomorphic( mi, mj );
    Sim sim = static_cast<Sim>( mode );
    mi = check_input( mi, sim );
    mj = check_input( mj, sim );
    SimFunc func = get_func( sim );
    RVecD out( mi.rows() );
    InterWrkr wrkr( mi, mj, out, func );
    RcppParallel::parallelFor( 0, mi.rows(), wrkr );
    return out;
}
