#include "innerp.hpp"

#if INNERP_IMPL != IMPL_CUDA

InnerP_func func_dot() {
    return []( Vec vi, Vec vj, ind i, ind j ) -> double {
        return vi.dot( vj );
    };
}

InnerP_func func_cosine() {
    return []( Vec vi, Vec vj, ind i, ind j ) -> double {
        return func_dot()( vi, vj, i, j ) / ( vi.norm() * vj.norm() );
    };
}

#endif
