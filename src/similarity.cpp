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

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends( RcppParallel )]]
// [[Rcpp::depends(RcppProgress)]]
#include "similarity.hpp"
#include "tools.hpp"
#include <functional>
#include <RcppParallel.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// using namespace Rcpp;

namespace impl {

    template<typename M1,typename M2>
    struct RwiseWorker : RcppParallel::Worker {
        const M1 a;
        const M2 b;
        const SimFunc func;
        RcppParallel::RVector<double> tgt;
        Progress &p;

        RwiseWorker( const M1 a, const M2 b, RVecD tgt, const SimFunc func, Progress &p ) :
            a( a ), b( b ), tgt( tgt ), func( func ), p( p ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( std::size_t i = begin; i < end; i++ ) {
                // if( p.check_abort() ) break;
                tgt[i] = func( a.row( i ), b.row( i ), i, i );
                // p.increment();
            }
        }
    };

    template<typename M1,typename M2>
    RVecD rwise_product( M1 a, M2 b, SimFunc func, bool quiet ) {
        isomorphic( a, b );
        RVecD c( a.rows() );
        Progress p( a.rows(), !quiet );
        RwiseWorker<M1,M2> wrkr( a, b, c, func, p );
        RcppParallel::parallelFor( 0, a.rows(), wrkr );
        return c;
    }

    template<typename M1,typename M2>
    struct InnerWorker : RcppParallel::Worker {
        const M1 a;
        const M2 b;
        const SimFunc func;
        RcppParallel::RMatrix<double> tgt;
        Progress &p;

        InnerWorker( const M1 a, const M2 b, RMatD tgt, const SimFunc func, Progress &p )
            : a( a ), b( b ), tgt( tgt ), func( func ), p( p ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( std::size_t i = begin; i < end; i++ ) {
                for( std::size_t j = 0; j < b.rows(); j++ ) {
                    // if( p.check_abort() ) break;
                    double v = func( a.row( i ), b.row( j ), i, j );
                    tgt( i, j ) = v;
                    p.increment();
                }
            }
        }
    };

    template<typename M1, typename M2>
    RMatD inner_product( M1 a, M2 b, SimFunc func, bool quiet ) {
        if( a.cols() != b.cols() ) Rcpp::stop( "Matrices must have equal column dimensions" );
        RMatD c( a.rows(), b.rows() );
        Progress p( a.rows() * b.rows(), !quiet );
        InnerWorker<M1,M2> wrkr( a, b, c, func, p );
        RcppParallel::parallelFor( 0, a.rows(), wrkr );
        return c;
    }

    template<typename M>
    struct SelfWorker : RcppParallel::Worker {
        const M src;
        const SimFunc func;
        const bool symm;
        RcppParallel::RMatrix<double> tgt;
        Progress &p;

        SelfWorker( const M a, RMatD tgt, const SimFunc func, const bool symm, Progress &p )
            : src( a ), tgt( tgt ), func( func ), symm( symm ), p( p ) {}

        void operator()( std::size_t begin, std::size_t end ) {
            for( std::size_t i = begin; i < end; i++ ) {
                for( std::size_t j = 0; j < ( symm ? ( i + 1 ) : src.rows() ); j++ ) {
                    // if( p.check_abort() ) break;
                    double v = func( src.row( i ), src.row( j ), i, j );
                    tgt( i, j ) = v;
                    if( symm ) tgt( j, i ) = v;
                    p.increment();
                }
            }
        }
    };

    template<typename M>
    RMatD self_product( M a, SimFunc func, bool symm, bool quiet ) {
        RMatD c( a.rows(), a.rows() );
        Progress p( a.rows() * a.rows(), !quiet );
        SelfWorker<M> wrkr( a, c, func, symm, p );
        RcppParallel::parallelFor( 0, a.rows(), wrkr );
        return c;
    }

    RowMat dirichlet_prior( RowMat m ) {
        double d = 1 / m.sum();
        return ( m.array() + d ).matrix().eval();
    }
}

//' Compute similarity or divergence/distance measures.
//'
//' This function will compute the designated similarity or distance function between corresponding
//' row (or column) vectors in the given matrices. The implementation assumes
//' observations are stored as rows, with features stored in columns, as per the usual R data frame
//' format; setting trans to TRUE will interpret the given matrices in the opposite direction.
//'
//' If no B matrix is given, the function will compute the measure between A and itself.
//'
//' If inner is set to FALSE, the function will compute the measure row-wise, returning a vector of
//' length equal to the number of rows in A and B where each entry i will correspond to the value
//' of the measure between A( i ) and B( i ).
//'
//' If B is given and inner is TRUE, A and B must have the same column dimension but any number of
//' rows. If B is given and inner is FALSE, A and B must be isomorphic.
//'
//' This function \emph{always} returns a dense matrix, as similarity and divergence measures
//' generally do not maintain the sparsity of the input. Because of this, applying this function
//' over a cooccurrence matrix containing feature vectors for every word in a lexicon is not
//' recommended and will most likely cause an out of memory error. The \emph{length} of the feature
//' vectors (i.e. the number of columns in A/B if trans is FALSE) on the other hand has no effect
//' over the size of the returned matrix, so the feature vectors themselves may cover the entire
//' lexical set. Do note however that running time will scalate linearly on the dimensionality of
//' the feature vectors (and quadratically on the number of observations).
//'
//' Available measure functions:
//' \itemize{
//'     \item{0: Additive cooccurrence retrieval. See \code{sim_additive}.}
//'     \item{1: Difference-weighted cooccurrence retrieval. See \code{sim_dweighted}.}
//'     \item{2: Cosine similarity. See \code{sim_cosine}.}
//'     \item{3: Kullback-Leibler divergence. See \code{div_kullback_leibler}.}
//'     \item{4: Jensen-Shannon divergence. See \code{div_jensen_shannon}.}
//'     \item{5: Bhattacharyya divergence. See \code{div_bhattacharyya}.}
//'     \item{6: Hellinger distance. See \code{dist_hellinger}.}
//' }
//'
//' Support for custom measures will be implemented pending API stabilization.
//'
//' @param A     A dense or sparse numeric matrix.
//' @param B     A dense or sparse numeric matrix with similar column dimension to A. NULL by
//'              default.
//' @param mode  The desired measure. See details.
//' @param inner Logical indicating whether to compute inner or row-wise measure. TRUE by default.
//'              Ignored if B is NULL.
//' @param trans Logical indicating whether the given matrices should be transposed before
//'              computation. Defaults to FALSE (i.e. compute measure over row vectors).
//'
//' @return If inner is TRUE and B is not null or if B is NULL, a dense matrix with as many rows as
//'         rows in A, and as many columns as there are rows in A if B is NULL or in B if it is not
//'         NULL. Each entry (i,j) will correspond to the value of the measure on A( i ) and B( j ),
//'         or if B is NULL, A( j ). If inner is FALSE, a dense vector of length equal to the
//'         number of rows in A and B where each entry i corresponds to the measured value between
//'         A( i ) and B( i ).
//'
//' @export
// [[Rcpp::export]]
Rcpp::RObject simdiv(
    Rcpp::RObject A,
    Rcpp::Nullable<Rcpp::RObject> B=R_NilValue,
    int mode=2, bool inner=true, bool trans=false, bool quiet=false, bool rm_zeros=false
) {
    Sim sim = static_cast<Sim>( mode );
    SimFunc func = get_func( sim );
    int a_type = A.sexp_type();
    Rcpp::RObject bm;
    Rcpp::RObject out = R_NilValue;
    if( B.isNull() ) {
        bm = A;
        if( a_type == REALSXP ) {
            RowMat a = asRowMajor<RowMat,Mat>( Rcpp::as<Mat>( A ), trans );
            a = rm_zeros ? impl::dirichlet_prior( a ) : a;
            out = impl::self_product( a, func, symmetric( sim ), quiet );
        } else if( a_type == S4SXP ) {
            RowSpMat a = asRowMajor<RowSpMat,SpMat>( Rcpp::as<SpMat>( A ), trans );
            if( rm_zeros ) Rcpp::warning( "Can't remove zeros from sparse matrix" );
            out = impl::self_product( a, func, symmetric( sim ), quiet );
        } else {
            Rcpp::stop( "Unsupported matrix type" );
        }
    } else {
        bm = Rcpp::as<Rcpp::RObject>( B.get() );
        int b_type = bm.sexp_type();
        if( a_type == REALSXP && b_type == REALSXP ) {      // Dense,Dense
            RowMat a = asRowMajor<RowMat,Mat>( Rcpp::as<Mat>( A ), trans );
            a = rm_zeros ? impl::dirichlet_prior( a ) : a;
            RowMat b = asRowMajor<RowMat,Mat>( Rcpp::as<Mat>( bm ), trans );
            b = rm_zeros ? impl::dirichlet_prior( b ) : b;
            if( inner ) {
                out = impl::inner_product( a, b, func, quiet );
            } else {
                out = impl::rwise_product( a, b, func, quiet );
            }
        } else if( a_type == S4SXP && b_type == REALSXP ) { // Sparse,Dense
            RowSpMat a = asRowMajor<RowSpMat,SpMat>( Rcpp::as<SpMat>( A ), trans );
            if( rm_zeros ) Rcpp::warning( "Can't remove zeros from sparse matrix" );
            RowMat b = asRowMajor<RowMat,Mat>( Rcpp::as<Mat>( bm ), trans );
            b = rm_zeros ? impl::dirichlet_prior( b ) : b;
            if( inner ) {
                out = impl::inner_product( a, b, func, quiet );
            } else {
                out = impl::rwise_product( a, b, func, quiet );
            }
        } else if( a_type == REALSXP && b_type == S4SXP ) { // Dense,Sparse
            RowMat a = asRowMajor<RowMat,Mat>( Rcpp::as<Mat>( A ), trans );
            a = rm_zeros ? impl::dirichlet_prior( a ) : a;
            RowSpMat b = asRowMajor<RowSpMat,SpMat>( Rcpp::as<SpMat>( bm ), trans );
            if( rm_zeros ) Rcpp::warning( "Can't remove zeros from sparse matrix" );
            if( inner ) {
                out = impl::inner_product( a, b, func, quiet );
            } else {
                out = impl::rwise_product( a, b, func, quiet );
            }
        } else if( a_type == S4SXP && b_type == S4SXP ) {   // Sparse,Sparse
            RowSpMat a = asRowMajor<RowSpMat,SpMat>( Rcpp::as<SpMat>( A ), trans );
            RowSpMat b = asRowMajor<RowSpMat,SpMat>( Rcpp::as<SpMat>( bm ), trans );
            if( rm_zeros ) Rcpp::warning( "Can't remove zeros from sparse matrix" );
            if( inner ) {
                out = impl::inner_product( a, b, func, quiet );
            } else {
                out = impl::rwise_product( a, b, func, quiet );
            }
        } else {
            Rcpp::stop( "Unsupported matrix format" );
        }
    }
    if( out != R_NilValue ) {
        if( A.hasAttribute( "Dimnames" ) ) {
            Rcpp::List anames = Rcpp::as<Rcpp::List>( A.attr( "Dimnames" ) );
            Rcpp::CharacterVector a_rows = Rcpp::as<Rcpp::CharacterVector>( anames[0] );
            if( !inner ) out.attr( "names" ) = a_rows;
            else {
                Rcpp::rownames( out ) = a_rows;
            }
        }
        if( bm.hasAttribute( "Dimnames" ) && inner ) {
            Rcpp::List bnames = Rcpp::as<Rcpp::List>( bm.attr( "Dimnames" ) );
            Rcpp::CharacterVector b_rows = Rcpp::as<Rcpp::CharacterVector>( bnames[0] );
            Rcpp::colnames( out ) = b_rows;
        }
    }
    return out;
}


