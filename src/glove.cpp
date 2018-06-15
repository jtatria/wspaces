#include "glove.hpp"

Glove::Glove(
    const SpMat& cooc, Mat& w, Mat& c, Vec& w_b, Vec& c_b, const Params params
) :
    m_params( params ), m_cooc( cooc ),
    m_w( w ),  m_c( c ), m_w_b( w_b ), m_c_b( c_b ),
    m_w_grad( w.rows(), w.cols() ),
    m_c_grad( c.rows(), c.cols() ),
    m_w_grad_b( w_b.size() ),
    m_c_grad_b( c_b.size() )
{
}

void Glove::epoch() {
    // TODO: port external C code to C++
}

Mat Glove:: w(){
    return m_w;
}

Mat Glove:: c(){
    return m_c;
}

Vec Glove:: w_bias(){
    return m_w_b;
}

Vec Glove:: c_bias(){
    return m_c_b;
}

Mat Glove:: wgrad(){
    return m_w_grad;
}

Mat Glove:: cgrad(){
    return m_c_grad;
}

Vec Glove:: wgrad_bias(){
    return m_w_grad_b;
}

Vec Glove:: cgrad_bias(){
    return m_w_grad_b;
}
