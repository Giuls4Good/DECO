// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// Correctly setup the build environment
// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
//Include deco.h AFTER namespace Rcpp!
#include "deco.h"

//' DECO Parallelized Algorithm (Pure C++)
//'
//' A description goes here
//'
//' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
//' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
//' @param p is the column dimension of X [equivalently, p is the number of regressor variables].
//' If not given, it is computed as the number of columns of X.
//' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals]
//'  If not given, it is computed as the number of rows of X.
//' @param m is the number of groups/blocks you wish to split X into, denoted X(i) for 1 <= i <= m
//' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
//' @param ncores determines the number of cores used on each machine to parallelize computation
//' @param r_1 is a tweaking parameter for making the inverse more robust (as we take inverse of XX + r_1*I)
//' @param r_2 is a tweaking parameter for making the inverse more robust (as we take inverse of X_MX_M + r_2*I)
//' @param intercept determines whether to include an intercept in the model or not
//' @param refinement determines whether to include the refinement step (Stage 3 of the algorithm)
//' @param parallel_lasso determines whether a parallel version of the Lasso coefficients should be used.
//' This parameter is ignored when \code{glmnet} is set to \code{FALSE} (see details).
//' @param glmnet determines whether glmnet function form glmnet R package should be used to compute the Lasso coefficients.
//' See details for further information. If set to \code{FALSE}, C++ implementation of coordinate descent algorithm is used.
//' @param precision determines the precision used in the coordinate descent algorithm. It is ignored when
//' \code{glmnet} is set to \code{TRUE}.
//' @param max_iter determines the maximum number of iterations used in the coordinate descent algorithm.
//' It is ignored when \code{glmnet} is set to \code{TRUE}.
//' @details This function is a pure C++ implementation of \code{DECO_LASSO_R} and \code{DECO_LASSO_MIX} functions.
//' Due to the fact that it is entirely written in C++ is generally way faster than its counterparts.
//'
//' Two functions can be used to compute Lasso coefficients: glmnet R function (\code{glmnet = TRUE})
//' and coordinate descent algorithm (\code{glmnet = FALSE}). glmnet R function is generally faster, but more memory is
//' required to pass the input argumentd from C++ to R and back. When \code{parallel_lasso = TRUE} an R parallelized
//' version of glmnet is used. Note however that for small datasets this could lead to slower run times, due to the
//' communication between C++ and R.
//'
//' Descent coordinate algorithm is always run on a single core.
//'
//' @return An estimate of the coefficients b.
//' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
//' @export
// [[Rcpp::export]]
arma::mat DECO_LASSO_C_PARALLEL(arma::vec& Y, arma::mat& X, int p, int n, int m, float lambda, float r_1, float r_2=0.01,
                        int ncores = 1, bool intercept = true, bool refinement = true, bool parallel_lasso = false,
                        bool glmnet = true, double precision = 0.0000001, int max_iter = 100000) {
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
  if(glmnet) {
    Environment RDeco("package:RDeco");
    if(parallel_lasso) {
      Function lassoCoefParallel = RDeco["lassoCoefParallel"];
      coefs = lassoRCoefParallel(X_new,Y_new,1.0,2*lambda,false,m,lassoCoefParallel);
    } else {
      Function lassoCoef = RDeco["lassoCoef"]; //Getting the function lassoCoef from package RDeco
      for(int i=0; i<m; i++) { //Can NOT parallelize this because it calls an R function!
        int startGroup = p_over_m*i;
        int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
        coefs.subvec(startGroup,endGroup) = lassoRCoef(X_new[i],Y_new,1.0,2*lambda,false,lassoCoef);
      }
    }
  } else {
    if((p-p_over_m*(m-1)) > n) { //Gradient descent does not work if p>n. p-p_over_m*(m-1) gives the greatest p
      //of all the Xi (note that the last Xi has necessarly the greates p of all)
      throw std::range_error("Gradient descent works only when all the m subsets have n<p. Use glmnet instead.");
    }
    //#pragma omp parallel for num_threads(ncores) -> CRASHES
    for(int i=0; i<m; i++) {
      int startGroup = p_over_m*i;
      int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
      coefs.subvec(startGroup,endGroup) = coordinateDescent_naive(X_new[i],Y_new,2*lambda,precision,max_iter);
    }
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
