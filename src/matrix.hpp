#include <RcppEigen.h>

#ifndef MATRIX_
#define MATRIX_ 1

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::Map<SpMat> MSpMat;

#endif // MATRIX_
