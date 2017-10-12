#ifndef TYPES_
#define TYPES_ 1

// [[Rcpp::plugins("cpp11")]]
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

typedef Rcpp::NumericMatrix RMatD;
typedef Rcpp::NumericVector RVecD;
typedef Rcpp::IntegerMatrix RMatI;
typedef Rcpp::IntegerVector RVecI;

typedef Eigen::Index ind;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::Map<SpMat> MSpMat;
typedef SpMat::InnerIterator SpInIt;

template <typename... args >
using F = std::function<double(args...)>;

#endif // TYPES_
