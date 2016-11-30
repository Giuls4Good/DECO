#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "deco.h"

//' Implement the refinement step in C++
//'
//' @param Y
//' @param X
//' @param coefs vector of estimated coefficients to be refined
//' @param intercept a parameter that if TRUE includes an intercept
//' @param r_2 a parameter that allows for matrix inversion to take place
//' @param n
//' @param p
//' @param lambda is the parameter in the lasso computations
//' @return coefs updated by the ridge regression
arma::vec refine(arma::vec &Y, arma::mat &X, arma::vec coefs, bool intercept, float r_2, int n, int p, float lambda){
  if(intercept){
    /*
    vec coef0 = mean(Y) - mean(X,0)*coefs;
    coefs.insert_rows(0,coef0);
     */
    mat column_of_ones = ones<mat>(X.n_rows);
    X = join_rows(column_of_ones, X);
  }

  //   STEP 4.1 Check if n<=#nonzero coefs and perform LASSO if so **#
  uvec M = find(abs(coefs) > 0.001);

  int number_of_non_zero_coefs = M.n_elem;
  //Subset X
  mat X_M = X.cols(M);

  if(number_of_non_zero_coefs >= n){
    Environment RDeco("package:RDeco");
    Function lassoCoef = RDeco["lassoCoef"]; //Getting the function lassoCoef from package mod4packcpp
    coefs.elem(M) = lassoRCoef(X_M,Y,1.0,2*lambda,false,lassoCoef);
    //coefs[M] = coef(glmnet(X[, M], Y, alpha = 1, nlambda = 1, lambda = c(2*lambda), intercept = FALSE))[-1];
    M = find(abs(coefs) > 0.001);
    X_M = X.cols(M);
  }

  if(number_of_non_zero_coefs == 0){
    printf("All coefficients estimated as 0.\n");
    //return(zeros<vec>size((p+1)));
    return(coefs);
  }
  //  STEP 4.2 Run Ridge regression on all non-zero coef vars     **#
  //Find the indicies of the coefficients that are non-zero

  //Apply Ridge Regression to give an updated and hopefully better estimate of the coefficient vector
  //coefs[M] = solve(t(X_M)%*%X_M + r_2*diag(length(M)))%*% t(X_M) %*% Y
  mat diagmat; diagmat.eye(M.size(),M.size());
  coefs.elem(M) = inv_sympd(X_M.t() * X_M + r_2*diagmat) * X_M.t() * Y;

  return(coefs);
}

//refine(1:3, matrix(1:12,3,4), c(1,2,-3,4), 1, r_2 = 0.001, 3, 4, 0.03)
