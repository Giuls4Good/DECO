#include <RcppArmadillo.h>
using namespace Rcpp;

//' Square root of a symmetric matrix
//'
//' @param M a symmetric matrix
//' @return A matrix C which is the square root of M (i.e. C*C=M)
//' @details To compute the square root of the matrix, sqrtmat_sympd function of Armadillo library is used.
//' @note This function is about 100 times faster than R function sqrtm (contained in expm package)
//' Note that no check is done to test if the matrix M is actually symmetric.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat squareRootSymmetric(arma::mat M) {
  arma::mat C = arma::sqrtmat_sympd(M);
  return(C);
}

//' Standardize a vector so that its mean is equal to 0
//'
//' @param V a vector
//' @return A vector with mean 0
//' @details The vector is standardized by subtracting each of its elements with the mean of V.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec standardizeVector(arma::vec V) {
  return(V-mean(V));
}

//' Standardize a matrix so that its mean is equal to 0
//'
//' @param M a matrix
//' @return A standardized matrix
//' @details The matrix is standardized by subtracting each column with colmeans (shitty explanation)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat standardizeMatrix(arma::mat M) {
  arma::mat mean_col = mean(M,0);
  for(size_t j=0; j<M.n_cols; j++) {
    M.col(j) -= mean_col(j);
  }
  return(M);
}

//' Transpose of a matrix
//'
//' @param M a matrix
//' @return Its transpose
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat tMatrix(arma::mat M) {
  return(M.t());
}

//' Multiply two matrices
//'
//' @param A a matrix
//' @param B a matrix whose size is compatible with the size of A
//' @return The matrix product A*B
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mulMatrices(arma::mat A, arma::mat B) {
  return(A*B);
}

//' Inverse of a matrix
//'
//' @param M a symmetric quadratic matrix
//' @return returns the inverse of the matrix
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat invSymmMatrix(arma::mat M) {
  return(inv_sympd(M));
}
