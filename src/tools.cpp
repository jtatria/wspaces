#include "tools.hpp"

Vec marginal( const SpMat &src, const Margin &margin, const F<double,double> fold ) {
    Vec out;
    ind (SpInIt::*index)() const;
    switch( margin ) {
        case Col : {
            out = Vec::Zero( src.cols() );
            index = &SpInIt::row;
            break;
        }
        case Row : {
            out = Vec::Zero( src.rows() );
            index = &SpInIt::col;
            break;
        }
        default : throw std::domain_error( "Invalid margin requested" );
    }

    for( int i = 0; i < src.outerSize(); i++ ) {
        for( SpInIt it( src, i ); it; ++it ) {
            ind tInd = (it.*index)();
            out[ tInd ] = fold( out[ tInd ], it.value() );
        }
    }
    return out;
}

Vec marg_sum( const SpMat &src, const Margin &margin ) {
    F<double,double> sum = [](double d1, double d2) { return d1 + d2; };
    return marginal( src, margin, sum );
}

Vec marg_prb( const SpMat &src, const Margin &margin ) {
    Vec cts = marg_sum( src, margin );
    return ( cts.array() / cts.sum() ).matrix();
}

Vec marg_avg( const SpMat &src, const Margin &margin ) {
    Vec cts = marg_sum( src, margin );
    return ( cts.array() / cts.size() ).matrix();
}

template<>
RowMat asRowMajor<RowMat,Mat>( Mat m, bool transpose ) {
    return RowMat( transpose ? m.transpose() : m );
}

template<>
RowSpMat asRowMajor<RowSpMat,SpMat>( SpMat m, bool transpose ) {
    return RowSpMat( transpose ? m.transpose() : m );
}

template<>
ColMat asColMajor<ColMat,Mat>( Mat m, bool transpose ) {
    return RowMat( transpose ? m.transpose() : m );
}

template<>
ColSpMat asColMajor<ColSpMat,SpMat>( SpMat m, bool transpose ) {
    return RowSpMat( transpose ? m.transpose() : m );
}
