#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat DECO_LASSO_C(arma::vec Y, arma::mat X, int p, int n, float lambda, float r,
                       int ncores = 1, bool intercept = true) {
  // STEP 0: INITIALIZATION OF VARIABLES
  arma::mat Y_stand, X_stand(X), mean_col(mean(X,0)), XX, X_new;
  arma::vec Y_new, coefs, coef0;
  arma::mat *Xi = new arma::mat[ncores];

  // STEP 1: INITIALIZATION

  // STEP 1.1 Mean Standardization
  //Standardize Y
  Y_stand  = Y-mean(Y);
  //Standardize X
  for(size_t j=0; j<X_stand.n_cols; j++) {
    X_stand.col(j) -= mean_col(j);
  }

  // STEP 1.2 Arbitrary Partitioning
  int p_over_m = p/ncores;
  int startGroup, endGroup;
  for(int i=0; i<ncores; i++) { //Do NOT parallelize, armadillo does not like it very much since you are accessing
                                //X_stand at each iteration -> Rstudio crashes
    startGroup = p_over_m*i;
    endGroup = (i!=(ncores-1)) ? p_over_m*(i+1)-1 : (p-1);
    Xi[i] = X_stand.submat(0, startGroup, n-1, endGroup);
  }

  // STEP 1.3 Distribution on machines
  //Only relevant in implementations that actually load the data into different machines

  // STEP 2: DECORRELATION

  // STEP 2.1 Compute X(i)'X(i) for each i  ###TO PARALLELALIZE ###

  //  STEP 2.2 Compute X'X from X(i)'X(i)  ###TO PARALLELALIZE###
  XX = X*X.t();

  // STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}
  XX.diag() += r;
  XX = sqrt(p)*arma::sqrtmat_sympd(inv_sympd(XX)); //Parallelizable (?)

  // STEP 2.4 Compute Y* and X*(i) for each i ###TO PARALLELALIZE###
  Y_new = XX*Y;
  X_new = XX*X;

  // STEP 3: ESTIMATION

  // STEP 3.1 LASSO for coefs on each partition i
  coefs.ones(p);
  delete[] Xi;

  // STEP 3.2 Put together all estimated coefs

  // STEP 3.3 Compute Intercept from coefs
  if(intercept) {
    coef0 = mean(Y) - mean(X,0)*coefs;
    coefs.insert_rows(0,coef0);
  }


  //  STEP 4: REFINEMENT

  return(coefs);
}
