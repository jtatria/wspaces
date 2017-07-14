// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

const int TF_BOOLEAN = 0;
const int TF_RAW     = 1;
const int TF_NORM    = 2;
const int TF_LOGNORM = 3;
const int TF_05NORM  = 4;
const int TF_KNORM   = 5;

const int IDF_UNARY  = 0;
const int IDF_PLAIN  = 1;
const int IDF_SMOOTH = 2;
const int IDF_MAX    = 3;
const int IDF_PROB   = 4;

Vec tf_bool(    const Vec &tfs, const double &L );
Vec tf_raw(     const Vec &tfs, const double &L );
Vec tf_norm(    const Vec &tfs, const double &L );
Vec tf_logNorm( const Vec &tfs, const double &L );
Vec tf_05norm(  const Vec &tfs, const double &L );
Vec tf_Knorm(   const Vec &tfs, const double &L, const double &K );

inline Vec tf( int &mode, const Vec &tfs, const double &L ) {
  Vec (*func)( const Vec&, const double& );
  switch( mode ) {
    case TF_BOOLEAN: func = tf_bool;    break;
    case TF_RAW:     func = tf_raw;     break;
    case TF_NORM:    func = tf_norm;    break;
    case TF_LOGNORM: func = tf_logNorm; break;
    case TF_05NORM:  func = tf_05norm;  break;
    default: Rcpp::stop( "Unknown mode for TF function" );
  }
  return (*func)( tfs, L );
}

Vec idf_unary(  const Vec &dfs, const double &D );
Vec idf_plain(  const Vec &dfs, const double &D );
Vec idf_smooth( const Vec &dfs, const double &D );
Vec idf_max(    const Vec &dfs, const double &D );
Vec idf_prob(   const Vec &dfs, const double &D );

inline Vec idf( const int &mode, const Vec &dfs, const double &D ) {
  Vec (*func)( const Vec&, const double &D );
  switch( mode ) {
    case 0: func = idf_unary; break;
    case 1: func = idf_plain; break;
    case 2: func = idf_smooth; break;
    case 3: func = idf_max; break;
    case 4: func = idf_prob; break;
    default: Rcpp::stop( "Unknown mode for IDF function" );
  }
  return (*func)( dfs, D );
}

inline Vec tf_bool( const Vec &tfs, const double &L ) {
  return tfs.cwiseSign();
}

inline Vec tf_raw( const Vec &tfs, const double &L ) {
  return tfs;
}

inline Vec tf_norm( const Vec &tfs, const double &L ) {
  return tfs / L;
}

inline Vec tf_logNorm( const Vec &tfs, const double &L ) {
  return ( ( tfs / L ).array() + 1 ).log().matrix();
}

inline Vec tf_05norm( const Vec &tfs, const double &L ) {
  return tf_Knorm( tfs, L, 0.5 );
}

inline Vec tf_Knorm( const Vec &tfs, const double &L, const double &K = 0.5 ) {
  return ( tfs / tfs.maxCoeff() ) * ( K + ( 1 - K ) );
}

inline Vec idf_unary( const Vec &dfs, const double &D ) {
  return dfs.cwiseSign();
}

inline Vec idf_plain( const Vec &dfs, const double &D ) {
  return ( dfs.array().inverse() * D ).log().matrix();
}

inline Vec idf_smooth( const Vec &dfs, const double &D ) {
  return ( ( dfs.array().inverse() * D ) + 1 ).log().matrix();
}

inline Vec idf_max( const Vec &dfs, const double &D ) {
  return ( dfs.array().inverse() * dfs.maxCoeff() ).log().matrix();
}

inline Vec idf_prob( const Vec &dfs, const double &D ) {
  return ( ( dfs.array().inverse() ) * ( ( dfs.array() ) * -1 + D ) ).log().matrix();
}

