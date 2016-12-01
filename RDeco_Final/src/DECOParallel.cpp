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
//' This implements the algorithm DECO which was introduced in "DECOrrelated feature space partitioning
//' for distributed sparse regression" by Wang, Dunson, and Leng (2016). It assumes that we take the lasso to be the penalized
//' regression scehme.
//'
//' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e.
//' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor.
//' @param p is the column dimension of X [equivalently, p is the number of regressor variables].
//' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals].
//' @param m is the number of groups/blocks you wish to split X into, denoted X(i) for 1 <= i <= m.
//' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm.
//' @param ncores determines the number of threads used on each machine to parallelize computation.
//' @param r_1 is a tweaking parameter for making the inverse more robust (as we take inverse of XX + r_1*I).
//' @param r_2 is a tweaking parameter for making the inverse more robust (as we take inverse of X_MX_M + r_2*I).
//' @param intercept determines whether to include an intercept in the model or not.
//' @param refinement determines whether to include the refinement step (Stage 3 of the algorithm).
//' @param parallel_glmnet determines whether a parallel version of the Lasso coefficients should be used.
//' This parameter is ignored when \code{glmnet} is set to \code{FALSE} (see details).
//' @param glmnet determines whether \code{glmnet} function form \code{glmnet} R package should be used to compute the Lasso coefficients.
//' See details for further information. If set to \code{FALSE}, C++ implementation of coordinate descent algorithm is used.
//' @param precision determines the precision used in the coordinate descent algorithm. It is ignored when
//' \code{glmnet} is set to \code{TRUE}.
//' @param max_iter determines the maximum number of iterations used in the coordinate descent algorithm.
//' It is ignored when \code{glmnet} is set to \code{TRUE}.
//' @details This function is a C++ implementation of \code{DECO_LASSO_R} and \code{DECO_LASSO_MIX} functions.
//' Due to the fact that it is entirely written in C++ it runs faster than the corresponding R implementations for sufficiently large matrices.
//'
//' Two functions can be used to compute Lasso coefficients: \code{glmnet} R function (\code{glmnet = TRUE}).
//' and coordinate descent algorithm (\code{glmnet = FALSE}). \code{glmnet} R function is generally faster, but more memory is
//' required to pass the input argumentd from C++ to R and back. When \code{parallel_glmnet = TRUE} an R parallelized
//' version of \code{glmnet} is used. Note however that for small datasets this could lead to slower run times, due to the
//' communication between C++ and R.
//'
//' Descent coordinate algorithm is always run in a parallel way (using \code{ncores} threads).
//'
//' @return An estimate of the coefficients b.
//' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
//' @export
// [[Rcpp::export]]
arma::mat DECO_LASSO_C_PARALLEL(arma::vec& Y, arma::mat& X, int p, int n, int m, float lambda, float r_1, float r_2=0.01,
                        int ncores = 1, bool intercept = true, bool refinement = true,
                        bool glmnet = true,bool parallel_glmnet = false, double precision = 0.0000001, int max_iter = 100000) {
  // STEP 0: INITIALIZATION OF VARIABLES
  /*
   * Y_stand will contain the standardized version of Y (mean(Y_stand) = 0). The original Y is still used to compute the intercept
   * X_stand will contain the standardized version of X (mean(X) = 0) and is initialazed to be equal to X.
   * mean_col contains the mean of each column of X
   * XX will contain the decorellation matrix
   * Y_new will contain the product between the docrellation matrix XX and Y_stand
   * coef0 will contain the value of the intercept
   * coefs will contain the final coefficients returned to the user
   * Xi is an array of matrices of length equal to m. Its elements are the subset of X_stand matrix.
   * XiXi_list is an array of matrices of length equal to m. Its elements are the symmetric matrixes obtained
   * by multiplying each elements of Xi with its transpose.
   * X_new is an array of matrices of length equal to m. Its element are the product between the decorrelation
   * matrix XX and each elements of Xi
   */
  arma::mat X_stand(X), mean_col(mean(X,0)), XX;
  arma::vec Y_stand, Y_new, coef0, coefs;
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
  /*
   * The big matrix X is partitioned in m matrices by columns.
   * Each submatrix contains n rows and floor(p/m) columns, except for the last submatrix which
   * may contain more columns if p is not divisible by m. Since m is usually small, the last submatrix
   * will have few columns more than the others (computation times are then equal for each thread).
   * EXAMPLE:
   * p = 10, m=4
   * Xi[0] = X[1:n,1:2]; Xi[1] = X[1:n,3:4]; Xi[2] = X[1:n,5:6]; Xi[3] = X[1:n,7:10]
   */
  int p_over_m = p/m;
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    int startGroup = p_over_m*i;
    int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
    Xi[i] = X_stand.submat(0, startGroup, n-1, endGroup);
  }

  // STEP 2: DECORRELATION

  // STEP 2.1 Compute X(i)X(i)' for each i
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    XiXi_list[i] = Xi[i]*Xi[i].t(); //Symmetric matrix
  }
  // STEP 2.2 Compute XX' from X(i)X(i)'
  // NOTE: Probably step 2.1 and 2.2 can be merged together in a better way.
  // reduce(XX:+) can not be used because + is an overloaded operator of the mat class and pragma
  // does not support it.
  XX.zeros(XiXi_list[0].n_rows,XiXi_list[0].n_cols);
  for(int i=0; i<m; i++) {
      XX = XX + XiXi_list[i];
  }
  delete[] XiXi_list; //Free memory!

  // STEP 2.3 Get (X'X = rI)^{-0.5}
  XX.diag() += r_1;
  XX = sqrt(p)*arma::sqrtmat_sympd(inv_sympd(XX)); //*_sympd functions are specific for symmetric matrices

  // STEP 2.4 Compute Y* and X*(i) for each i
  Y_new = XX*Y_stand;
  #pragma omp parallel for num_threads(ncores)
  for(int i=0; i<m; i++) {
    X_new[i] = XX*Xi[i];
  }
  delete[] Xi; //free memory

  // STEP 3: ESTIMATION

  // STEP 3.1 LASSO for coefs on each partition i
  /*
   * If glmnet = true, parallel_glmnet = true -> R glmnet function in a parallel way is used
   * If glmnet = true, parallel_glmnet = false -> R glmnet function on a single core is used
   * If glmnet = false -> C++ gradient descent algorithm is used in a parallel way
   */
  if(glmnet) {
    Environment RDeco("package:RDeco"); //Get the environment defined by our package RDeco
    if(parallel_glmnet) {
      Function lassoCoefParallel = RDeco["lassoCoefParallel"]; //Get the R function that computes the Lasso
      //coefficients in a parallel way
      coefs = lassoRCoefParallel(X_new,Y_new,1.0,2*lambda,false,m,lassoCoefParallel); //Call the C++ function
      //that will call the R function.
    } else {
      Function lassoCoef = RDeco["lassoCoef"]; //Getting the function lassoCoef from package RDeco
      for(int i=0; i<m; i++) { //Can NOT parallelize this because it calls an R function!
        int startGroup = p_over_m*i;
        int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
        coefs.subvec(startGroup,endGroup) = lassoRCoef(X_new[i],Y_new,1.0,2*lambda,false,lassoCoef); //Call
        //the C++ function that calls the R function.
      }
    }
  } else {
    if((p-p_over_m*(m-1)) > n) { //Gradient descent does not work if p>n. p-p_over_m*(m-1) gives the greatest p
      //of all the Xi (note that the last Xi has necessarly the greates p of all)
      throw std::range_error("Gradient descent works only when all the m subsets have n<p. Use glmnet instead.");
    }
    #pragma omp parallel for num_threads(ncores)
    for(int i=0; i<m; i++) {
      int startGroup = p_over_m*i;
      int endGroup = (i!=(m-1)) ? p_over_m*(i+1)-1 : (p-1); //shortcut for if
      coefs.subvec(startGroup,endGroup) = coordinateDescent_naive(X_new[i],Y_new,lambda,precision,max_iter);
    }
  }
  delete[] X_new; //Free memory!

  // STEP 3.3 Compute Intercept from coefs
  if(intercept) {
    coef0 = mean(Y) - mean(X,0)*coefs;
    coefs.insert_rows(0,coef0);
  }

  //  STEP 4: REFINEMENT
  if(refinement) {
    coefs = refine(Y_stand, X_stand, coefs, intercept, r_2, n, p, lambda);
  }
  else if(intercept) {

  }
  return(coefs);
}
