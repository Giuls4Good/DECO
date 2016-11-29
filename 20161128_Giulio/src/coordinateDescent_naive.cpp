// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>     //e.g., ceil()

using namespace Rcpp; using namespace arma;

//CONTAINS: update_naive, coordinateDescent_naive


//' TITLE OF THIS FUNCTION
//'
//' @param beta gives the arma::vec vector of beta-estimates at the current iteration of Coordinate descent
//' @param X gives the arma::mat matrix of regressors in this partition of all regressors (there are m partitions
//'        in total, and we run coordinate descent on each partition)
//' @param Y gives the arma::vec vector of observations that we project X onto
//' @param lambda gives the double giving our lambda coefficient
//' @param n gives an integer corresponding to the number of rows (observations) of X (Y)
//' @param p gives the integer corresponding to the number of columns/regressors contained in X
//' @export
//' @return an updated version of beta
//'         NEEDS CHANGING:
//'         - we should compute Xi'Xi once for each i, and then supply it as another parameter?
//'         - we should be able to just return a pointer to a (changed) beta (more efficient than local copies!)
// [[Rcpp::export]]
arma::vec update_naive(arma::vec &beta, arma::mat &X, arma::vec &Y, double lambda, int n, int p){
  //THIS FUNCTION DEPENDA ON ARMADILLO!
  //STEP 1: declare the quantities needed
  double Soft_Thresholding_argument; //will hold the argument inside S(., lambda) (thresholding function)
  vec Xi;   //later on, this vector will hold X_i

  //STEP 2: LOOP to perform one gradient descent step
  for(int i = 0; i<p;i++){

    //Subsetting operations to obtain X and beta with 0s at the right places
    Xi = (X).col(i);    //initialize Xi with the first column vector of X
    //NOTE: Armadillo makes a local copy, so there is no danger in modifying Xi!
    //      I.e., it will not affect X if we modify Xi

    (X).col(i) = zeros<vec>(n);   //The first column is replaced by zeros!
    beta(i) = 0;                //Similarly, the first entry of beta is put to zero
    //This is done such that X%*%beta = X_-i%*%beta_-i!

    //compute X_-i%*%beta_-i
    Soft_Thresholding_argument = dot(Xi, (Y - (X*(beta)))) * (1.0/((double)n)); /* Xi_Squared_inv;  */
    //IMPLEMENT C ROUTINE ?

    //do the cases and update beta(i) using the thresholding

    if(Soft_Thresholding_argument > 0 && lambda < fabs(Soft_Thresholding_argument)){
      beta(i) = Soft_Thresholding_argument - lambda;
    }else if(Soft_Thresholding_argument < 0 && lambda < fabs(Soft_Thresholding_argument)){
      beta(i) = Soft_Thresholding_argument + lambda;
    }else{
      beta(i) = 0;
    }

    //put back Xi at the right place
    (X).col(i) = Xi;
  }

  //return the beta-update
  return(beta);

}


//' TITLE OF THIS FUNCTION
//'
//' @param X gives the arma::mat matrix of regressors in this partition of all regressors (there are m partitions
//'        in total, and we run coordinate descent on each partition)
//' @param Y gives the arma::vec vector of observations that we project X onto
//' @param lambda gives the double giving our lambda coefficient
//' @param precision gives the convergence criterion (how close two subsequent iterations should be before termination)
//' @param max_iter gives the maximum number of iterations in the inner update loop of coordinate descent
//' @export
//' @return an updated version of beta
//'         NEEDS CHANGING:
//'         - we should compute Xi'Xi once for each i, and then supply it as another parameter?
//'         - we should be able to just return a pointer to a (changed) beta (more efficient than local copies!)
// [[Rcpp::export]]
arma::vec coordinateDescent_naive(arma::mat X, arma::vec Y, double lambda, double precision, int max_iter){
  //THIS FUNCTION DEPENDA ON ARMADILLO!
  //STEP 1: INITIALIZATION
  unsigned int p = X.n_cols, n = X.n_rows; //get p (n) as the number of rows (columns) of the regression matrix
  vec beta, beta_update; //declare beta

  if(n != Y.n_elem){  //check if n = Y.n_elem and throw error if not
    throw std::range_error("row number of regressor matrix unequal to number of elements in Y for some partition!");
  }

  if(n>p){ //initialize using OLS or at random if p>=n
    beta = solve(X,Y);  // armadillo function x=solve(A,b) solves the system of linear equations A%*%x = b
    // this is faster than using .i() to get the inverse function of X'X
    //IMPLEMENT C ROUTINE
    //beta = randn<vec>(p);

  }else{
    beta = randn<vec>(p);   //initialize the coefficients as standard normaXi_Squared_inv, l random variables
    //if the OLS-equation A%*%x = b cannot be solved due to singularity of A=X'X.
  }

  //possibly: compute the inverse variances (should be passed on to update) - unclear if beneficial! Probably only for
  //very large n and p (due to having to store all of this in the memory, it might be quicker to recompute)

  //STEP 2: LOOP UNTIL CONVERGENCE ACHIEVED
  //loop until convergence criterion achieved
  bool criterion = false;
  int count = 0;
  while(!criterion && count<=max_iter){
    Rcpp::checkUserInterrupt();                                         //check if User pressed 'Stop'
    beta_update = update_naive(beta, X, Y, lambda, n, p);                 //get updated beta
    criterion =  ( (max(abs(beta_update - beta))) < precision );       //check if your convergence criterion is satisfied
    count++;                                                        //update the loop counter
    beta = beta_update;  //Update the beta
  }

  //return the beta-vector
  return(beta);
}
