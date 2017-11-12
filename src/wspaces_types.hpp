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
typedef Eigen::Map<Mat>   MMat;
typedef Eigen::Map<Vec>   MVec;
typedef SpMat::InnerIterator SpInIt;
typedef Mat::Scalar scalar;

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> ColMat;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMat;

typedef Eigen::SparseMatrix<double,Eigen::ColMajor,ind> ColSpMat;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor,ind> RowSpMat;

template <typename... args >
using F = std::function<scalar(args...)>;

#endif // TYPES_
