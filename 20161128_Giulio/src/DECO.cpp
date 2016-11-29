// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "deco.h" //Put it here, after all the previous stuff!

arma::vec lassoRCoef(arma::mat& X, arma::vec& Y, double nlambda,
                         double lambda, bool intercept, Function& f) { //The Function f is lassoCoef.
  NumericVector res = f(X, Y, nlambda, lambda, intercept);
  arma::vec res_arm(res.begin(), res.size(), false); //converting numericvector to armadillo by reference (not copying)
  return res_arm;
}

arma::vec lassoRCoefParallel(arma::mat* X, arma::vec& Y, double nlambda,
                        double lambda, bool intercept, int m, Function& f) {
  //Convert data Xi to a list so that it can be feed into R
  Rcpp::List data(m);
  for(int i=0; i<m; i++){
    data[i] = X[i];
  }
  NumericVector res = f(data, Y, nlambda, lambda, intercept, m);
  arma::vec res_arm(res.begin(), res.size(), false); //converting numericvector to armadillo by reference (not copying)
  return res_arm;
}

//' DECO Parallelized Algorithm (Pure C)
//'
//' This function is deprecated. Use \code{DECO_LASSO_C_PARALLEL} function.
//'
//' @details This function is equivalent to \code{DECO_LASSO_C_PARALLEL} function when fixing \code{m=1, ncores=1}.
//'
//' @export
// [[Rcpp::export]]
arma::mat DECO_LASSO_C(arma::vec Y, arma::mat X, int p, int n, float lambda, float r,
                       int ncores = 1, bool intercept = true) {
  throw std::invalid_argument("This function is defunct. Use DECO_LASSO_C_PARALLEL insted.");
  /*
  if(ncores > 1) {
    stop("This function does not support parallel computing. Check DECO_LASSO_C_PARALLEL function.");
  }
  // STEP 0: INITIALIZATION OF VARIABLES
  arma::mat Y_stand, X_stand(X), mean_col(mean(X,0)), XX, X_new;
  arma::vec Y_new, coefs, coef0;
  // STEP 1: INITIALIZATION

  // STEP 1.1 Mean Standardization
  //Standardize Y
  Y_stand  = Y-mean(Y);
  //Standardize X
  for(size_t j=0; j<X_stand.n_cols; j++) {
    X_stand.col(j) -= mean_col(j);
  }

  // STEP 1.2 Arbitrary Partitioning ###TO PARALLELIZE###

  // STEP 1.3 Distribution on machines
  //Only relevant in implementations that actually load the data into different machines

  // STEP 2: DECORRELATION

  // STEP 2.1 Compute X(i)'X(i) for each i  ###TO PARALLELALIZE ###

  //  STEP 2.2 Compute X'X from X(i)'X(i)  ###TO PARALLELALIZE###
  XX = X_stand*X_stand.t();
  // STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}
  XX.diag() += r;
  XX = sqrt(p)*arma::sqrtmat_sympd(inv_sympd(XX)); //Parallelizable (?)
  // STEP 2.4 Compute Y* and X*(i) for each i ###TO PARALLELALIZE###
  Y_new = XX*Y_stand;
  X_new = XX*X_stand;

  // STEP 3: ESTIMATION

  // STEP 3.1 LASSO for coefs on each partition i ###TO PARALLELALIZE ###
  Environment RDeco("package:RDeco");
  Function lassoCoef = RDeco["lassoCoef"]; //Getting the function lassoCoef from package RDeco
  coefs = lassoRCoef(X_new,Y_new,1.0,2*lambda,false,lassoCoef); //Calling C++ function that calls lassoCoef R function
  //coefs = lassoCoef(X_new,Y_new,1.0,2*lambda,false); ERROR
  // STEP 3.2 Put together all estimated coefs

  // STEP 3.3 Compute Intercept from coefs
  if(intercept) {
    coef0 = mean(Y) - mean(X,0)*coefs;
    coefs.insert_rows(0,coef0);
  }


  //  STEP 4: REFINEMENT
  return(coefs);
   */
}
