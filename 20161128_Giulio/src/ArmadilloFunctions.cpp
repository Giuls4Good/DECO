// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Square root of a symmetric matrix
//'
//' @param M a symmetric matrix
//' @return A matrix C which is the square root matrix of M (i.e. C*C=M)
//' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
//' @details To compute the square root of the matrix, \code{sqrtmat_sympd} function of Armadillo library is used.
//' @note This function is about 100 times faster than R function \code{sqrtm} (contained in expm package).
//'
//' Note that no check is done to test if the matrix M is actually symmetric.
//' @examples
//' require(expm)
//' require(rbenchmark)
//' A <- matrix(rnorm(10000,mean=10,sd=5),nrow=100)
//' A_symm <- A%*%t(A)
//' benchmark(sqrtm(A_symm),squareRootSymmetric(A_symm), order='relative')
//'
//' @export
// [[Rcpp::export]]
arma::mat squareRootSymmetric(arma::mat M) {
  arma::mat C = arma::sqrtmat_sympd(M);
  return(C);
}

//' Standardize a vector
//'
//' @param V a vector
//' @return Returns a vector with mean equal to 0.
//' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
//' @details To compute the standardized vector, each of its entries is subtracted with the vector's mean.
//' @note In general, this function does not return a vector with variance equal to 1.
//'
//' This function is about 2x slower than directly computing \code{v-mean(v)}, but it is faster than doing
//' \code{scale(v,scale=FALSE)}
//' @examples
//' require(rbenchmark)
//' v <- 1:5000000
//' benchmark({v-mean(v)},standardizeVector(v),scale(v,scale=FALSE), order='relative')
//'
//' @export
// [[Rcpp::export]]
arma::vec standardizeVector(arma::vec V) {
  return(V-mean(V));
}

//' Standardize a matrix so that its mean is equal to 0
//'
//' @param M a matrix
//' @return A matrix with mean 0
//' @details The matrix is standardized by subtracting each column with the mean of that column.
//' @note In general, this function does not return a matrix with variance equal to 1.
//'
//' This function is about 2.5x times faster than doing \code{scale(M,scale=FALSE)}.
//' @examples
//' require(rbenchmark)
//' M <- matrix(rnorm(1000*5000,10,5), nrow=1000)
//' benchmark(scale(M,scale=FALSE),standardizeMatrix(M), order='relative')
//'
//' @export
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
//' @details To compute the transpose matrix, Armadillo library is used.
//' @note This function is about 4 times slower than R function \code{t}.
//' @examples
//' require(rbenchmark)
//' benchmark(t(M),tMatrix(M),order='relative')
//'
//' @export
// [[Rcpp::export]]
arma::mat tMatrix(arma::mat M) {
  return(M.t());
}

//' Multiply two matrices
//'
//' @param A a matrix
//' @param B a matrix whose size is compatible with the size of A
//' @return The matrix product A*B.
//' @details To compute such product, Armadillo library is used.
//' @note This function takes about the same speed as R matrix multiplication.
//' @examples
//' require(rbenchmark)
//' A <- matrix(rnorm(1000*500,10,5), nrow=1000, ncol=500)
//' B <- matrix(rnorm(1000*500,10,5), nrow=500, ncol=1000)
//' benchmark(A%*%B,mulMatrices(A,B),order='relative')
//'
//' @export
// [[Rcpp::export]]
arma::mat mulMatrices(arma::mat A, arma::mat B) {
  return(A*B);
}

//' Inverse of a matrix
//'
//' @param M a symmetric quadratic matrix
//' @return The inverse of the matrix
//' @details To compute the square root of the matrix, \code{inv_sympd} function of Armadillo library is used.
//' @note This function is about 2.5x times faster than R function \code{solve}.
//'
//' Note that no check is done to test if the matrix M is actually symmetric.
//' @examples
//' require(rbenchmark)
//' M <- matrix(rnorm(1000^2,10,5), nrow=1000)
//' M_symm <- M%*%t(M)
//' benchmark(solve(M_symm),invSymmMatrix(M_symm),order='relative')
//'
//' @export
// [[Rcpp::export]]
arma::mat invSymmMatrix(arma::mat M) {
  return(inv_sympd(M));
}
