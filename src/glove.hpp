#ifndef _GLOVE_HPP
#define _GLOVE_HPP 1

#include "wspaces_types.hpp"

class Glove {

public:
    struct Params {
        const scalar m_alpha;
        const scalar m_rate;
        const scalar m_lambda;
        const int m_max_x;

        Params( const scalar alpha, const scalar rate, const scalar lambda, const int max_x )
        : m_alpha( alpha ), m_rate( rate ), m_lambda( lambda ), m_max_x( max_x ) {}
    };

    Glove( const SpMat& cooc, Mat& W_vs, Mat& C_vs, Vec& W_b, Vec& C_b, const Params params );

    void epoch();
    Mat w();
    Mat c();
    Vec w_bias();
    Vec c_bias();
    Mat wgrad();
    Mat cgrad();
    Vec wgrad_bias();
    Vec cgrad_bias();

private:
    const SpMat m_cooc;
    const Params m_params;

    Mat m_w, m_c, m_w_grad, m_c_grad;
    Vec m_w_b, m_c_b, m_w_grad_b, m_c_grad_b;
};

#endif
