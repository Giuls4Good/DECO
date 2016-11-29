// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]
// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
//Include deco.h AFTER namespace Rcpp!
#include "deco.h"

// [[Rcpp::export]]
arma::mat DECO_LASSO_C_PARALLEL(arma::vec Y, arma::mat X, int p, int n, int m, float lambda, float r_1,
                        int ncores = 1, bool intercept = true, bool refinement = true) {
  // STEP 0: INITIALIZATION OF VARIABLES
  arma::mat Y_stand, X_stand(X), mean_col(mean(X,0)), XX;
  arma::vec Y_new, coef0, coefs;
  arma::mat *Xi = new arma::mat[m]; //REMEMBER to delete it after use!!
  arma::mat *XiXi_list = new arma::mat[m]; //REMEMBER to delete it after use!!
  arma::mat *X_new = new arma::mat[m]; //REMEMBER to delete it after use!!
  coefs.ones(p); //Initialize dimension of p and fill with ones

  // STEP 1: INITIALIZATION

  // STEP 1.1 Mean Standardization
  //Standardize Y
  Y_stand  = Y-mean(Y);
  //Standardize X
  for(size_t j=0; j<X_stand.n_cols; j++) {
    X_stand.col(j) -= mean_col(j);
  }

  // STEP 1.2 Arbitrary Partitioning
  int p_over_m = p/m;
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    int startGroup = p_over_m*i;
    int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
    Xi[i] = X_stand.submat(0, startGroup, n-1, endGroup);
  }

  // STEP 1.3 Distribution on machines
  //Only relevant in implementations that actually load the data into different machines

  // STEP 2: DECORRELATION

  // STEP 2.1 Compute X(i)X(i)' for each i
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    XiXi_list[i] = Xi[i]*Xi[i].t();
  }
  // STEP 2.2 Compute XX' from X(i)X(i)'
  // Probably step 2.1 and 2.2 can be merged together in a better way
  XX.zeros(XiXi_list[0].n_rows,XiXi_list[0].n_cols);
  for(int i=0; i<m; i++) {
      XX = XX + XiXi_list[i];
  }
  delete[] XiXi_list;

  // STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}
  XX.diag() += r_1;
  XX = sqrt(p)*arma::sqrtmat_sympd(inv_sympd(XX)); //parallelizable??

  // STEP 2.4 Compute Y* and X*(i) for each i
  Y_new = XX*Y_stand;
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    X_new[i] = XX*Xi[i];
  }
  delete[] Xi; //free memory

  // STEP 3: ESTIMATION

  // STEP 3.1 LASSO for coefs on each partition i
  Environment RDeco("package:RDeco");
  Function lassoCoef = RDeco["lassoCoef"]; //Getting the function lassoCoef from package RDeco
  for(int i=0; i<m; i++) { //Can NOT parallelize this because it calls an R function!
    int startGroup = p_over_m*i;
    int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
    coefs.subvec(startGroup,endGroup) = lassoRCoef(X_new[i],Y_new,1.0,2*lambda,false,lassoCoef);
  }
  delete[] X_new;

  // STEP 3.3 Compute Intercept from coefs
  if(intercept) {
    coef0 = mean(Y) - mean(X,0)*coefs;
    coefs.insert_rows(0,coef0);
  }

  //  STEP 4: REFINEMENT
  if(refinement) {
    Rcpp::warning("Sorry, refinement step is still not implemented");
  }
  return(coefs);
}
