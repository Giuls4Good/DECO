#include <RcppArmadillo.h>
#include <math.h>     //e.g., ceil()

//#include <armadillo.h>
//using namespace Rcpp;

using namespace Rcpp; using namespace arma;
//using namespace std;
//using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//Useful commands:
//- Rprinft("...") prints ... to the R console


//a function to investigate how the soft thresholding computation is performed and which results it has
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec SThres(arma::vec beta, arma::mat X, arma::vec Y,double lambda, int n, int p){
  int i=0;

  vec Xi = X.col(i);
  X.col(i) = zeros<vec>(n);
  beta(i) = 0;

  double Soft_Thresholding_argument = dot(Xi, (Y - (X*beta))) * (1.0/((double) n));

  if(Soft_Thresholding_argument > 0 && lambda < fabs(Soft_Thresholding_argument)){
    beta(i) = Soft_Thresholding_argument - lambda;
  }else if(Soft_Thresholding_argument < 0 && lambda < fabs(Soft_Thresholding_argument)){
    beta(i) = Soft_Thresholding_argument + lambda;
  }else{
    beta(i) = 0;
  }


  vec res = zeros<vec>(2);
  res(0) = Soft_Thresholding_argument;
  res(1) = beta(i);
  X.col(i) = Xi;
  return(res);
}

//a more efficient way of updating
//' @param beta gives the arma::vec vector of beta-estimates at the current iteration of Coordinate descent
//' @param X gives the arma::mat matrix of regressors in this partition of all regressors (there are m partitions
//'        in total, and we run coordinate descent on each partition)
//' @param Y gives the arma::vec vector of observations that we project X onto
//' @param lambda gives the double giving our lambda coefficient
//' @param n gives an integer corresponding to the number of rows (observations) of X (Y)
//' @param p gives the integer corresponding to the number of columns/regressors contained in X
//' @param k gives the number of coefficients currently in the model (i.e., nonzero)
//' @param n_inv = 1/n. Passed on because it is needed frequently and one may avoid division this way
//' @param XXCP Cross-products between the columns of X that are active. Cross-product between X_i and X_j is stored in position
//'        (min(i,j), max(i,j)) (so we have an upper triangular matrix). Inactive regressors have zeros
//' @param XYCP Cross-products between X_i and Y at position i
//' @param inModelIndices gives a vector with p entries, of which the first k have meaning as the indices of active regressors
//' @param inModelBool gives a vector with p entries which are either 1 (true) or 0 (false) and indicate if the coef is in the model
//'         (unsigned char is used as it is the smallest data type available in armadillo)
//' @param gradients gives the vector of the p gradient values of the coefficients
//' @return an updated version of beta
//' @details based on the algorithm outlined in 'Regularization Paths for GLM via Coordinate Descent', Friedman, Hastie, Tibshirani, 2010.
//'          In particular, this algorithm uses the covariance-updating scheme laid out in section 2.2 of that paper
//'          NOTE: NEEDS CHANGING -> inModelIndices is not needed and probably slows down the code because of the deletion! (O(k))
//'                               -> give n_inv as parameter into the function
//'                               -> change (i,j) to .at(i,j) in order to avoid the checks and make the operation quicker
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::vec update_cov(arma::vec &beta, arma::mat &X, const arma::vec &Y, double lambda, unsigned int n, unsigned int p, unsigned int k,
                     double n_inv, arma::mat &XXCP, const arma::vec &XYCP, arma::Col<unsigned int> &inModelIndices,
                     arma::Col<unsigned char> &inModelBool, arma::vec &gradients){

  //STEP 1: declare quantities needed
  bool inModel_start, inModel_end, beta_i_old, helper;    //keep track if beta_i was in the model at the beginning and the end of the inner loop
                                              //and define a helper for finding the index i in inModelIndices
  double ThresholdArgument;
  unsigned int ind;


  //STEP 2: LOOP over all regressors
  for(unsigned int i=0;i<p;++i){

    //STEP 2.1 check if beta_i currently in model and get its current value as beta_i_old
    inModel_start = inModelBool(i);
    beta_i_old = beta(i);

    //STEP 2.2 perform Thresholding on the gradient and see if the coefficient is in the model at the end of it
    ThresholdArgument = gradients(i);

    if(ThresholdArgument > lambda){                        // ...   ...  [-lambda  ... 0 ... lambda]   ...  x  ...
      beta(i) = ThresholdArgument - lambda;
      inModel_end=1; //indicate that regressor i in model
    }else if(ThresholdArgument < (-lambda)){               // ... x ...  [-lambda  ... 0 ... lambda]   ...     ...
      beta(i) = ThresholdArgument + lambda;
      inModel_end=1; //indicate that regressor i in model
    }else{                                                 // ...   ...  [-lambda  ... x ... lambda]   ...     ...
      beta(i) = 0;
      inModel_end=0; //indicate that regressor i NOT in model
    }


    //STEP 2.3 check if the coefficient was put into the model this iteration for the first time.
    //         If so, fill XXCP with the relevant cross products and update inModelBool, inModelIndices, k

    if(inModel_start<inModel_end){

      //2.3.1 fill in the missing cross-products if beta_i was added in STEP 2.2
      for(unsigned int j=0; j<k; ++j){
        ind = inModelIndices(j);            //get the next index that belongs to an active regressor
        if( abs(XXCP(i,ind))>0.000001 ){    //update the cross-product matrix only if that wasn't done already
          XXCP(std::min(i,ind), std::max(i,ind)) = dot(X.col(i), X.col(ind))*n_inv;   //update the matrix entry
        }
      }

      //2.3.2 update inModelBool, inModelIndices, k
      inModelBool(i) = 1;     //so the i-th coefficient is part of the model now
      inModelIndices(k) = i;  //so the coefficient i is now one of the model coefficients
      k = k+1;                //we have another coefficient as part of the model!

    }

    //STEP 2.4 check if the coefficient left the model this iteration. If so, update inModelBool, inModelIndices, k

    if(inModel_start>inModel_end){
      inModelBool(i) = 0;     //so the i-th coefficient is part of the model now

      //2.4.1 find the position of i in inModelIndices
      helper=false; ind=-1;              //put ind equal to -1 for debugging purposes
      for(unsigned int j=0; (!helper && j<k); j++){

        if(inModelIndices(j) == i){
            helper=true;                //make sure we will exit the loop
            ind = j;                    //make sure we get the index of the desired value
        }

      }

      //2.4.2 update inModelBool, inModelIndices, k
      inModelIndices(ind) = inModelIndices(k-1);  //overwrite the position of the index leaving the model
      inModelIndices(k-1) = -1;                   //for debugging purposes
      k = k-1;                                  //we lost a coefficient that previously was part of the model!
    }


    //STEP 2.5 check if beta_i has changed (i.e.,  inModel_end != inModel_start || inModel_end == 1 ),
    //         where the first condition means that the role of beta_i has changed and the second condition means
    //         that the coefficient is (still) in the model.
    //         If the conditions are true, update the gradients

    if(inModel_start != inModel_end || inModel_end == 1){  //we need both conditions. Suppose you only use the first one: You miss the case
                                                          //where beta_i changes, but stays in the model. Suppose you only use the second one:
                                                         //You miss the case where beta_i left the model
      double diff = beta_i_old - beta(i);    //difference between beta_i before and after STEP 2.2
      double diff_x_ninv = diff*n_inv;       //difference scaled by 1/n

      //We need to update ALL gradients by adding <x_i, x_j> * (beta_i_old - beta_i_new) = XXCP(min(i,j), max(i,j))*diff
      for(unsigned int j=0; j<i; j++){       //Loop over ALL gradients EXCEPT for the i-th one (the one changed last)
                                    //Loop part 1: where j<i (because of how we store in XXCP, this split loop is beneficial)

        gradients(j) = gradients(j) + XXCP(j,i)*diff_x_ninv;  //here: j<i, so we do not need min/max computation

      }
      gradients(i) = gradients(i) + diff; //by standardization of the row vectors of X, this is possible
      for(unsigned int j=i;j<p; j++){     //Loop over ALL gradients EXCEPT for the i-th one (the one changed last)
                                   //Loop part 2: where j>i (because of how we store in XXCP, this split loop is beneficial)

        gradients(j) = gradients(j) + XXCP(i,j)*diff_x_ninv;  //here: j>i, so we do not need min/max computation

      }

    }


  }//end of for loop

  //return the parameter vector
  return(beta);
}

//THIS FUNCTION DEPENDA ON ARMADILLO!
//' @param X gives the arma::mat matrix of regressors in this partition of all regressors (there are m partitions
//'        in total, and we run coordinate descent on each partition)
//' @param Y gives the arma::vec vector of observations that we project X onto
//' @param lambda gives the double giving our lambda coefficient
//' @param precision gives the convergence criterion (how close two subsequent iterations should be before termination)
//' @param max_iter gives the maximum number of iterations in the inner update loop of coordinate descent
//' @param XXCP Cross-products between the columns of X that are active, see update_cov. Only passed on for memory reasons.
//' @param XYCP Cross-products between X_i and Y at position i, see update_cov. Only passed on for memory reasons.
//' @param inModelIndices gives a vector with p index entries, see update_cov. Only passed on for memory reasons.
//' @param inModelBool gives a vector with p entries which are either 1 (true) or 0 (false), see update_cov. Only passed on for memory reasons.
//' @param gradients gives the vector of the p gradient values of the coefficients, see update_cov. Only passed on for memory reasons.
//' @export
//' @return beta
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec coordinateDescent_cov(arma::mat &X, arma::vec &Y, double lambda, double precision, unsigned int max_iter){

  //STEP 1: INITIALIZATION

  //STEP 1.1: Create objects and variables needed later on
  unsigned int p = X.n_cols, n = X.n_rows, k, ind1, ind2; //get p (n) as the number of rows (columns) of the regression matrix
  vec beta, beta_update; //declare beta
  double n_inv = 1.0/((double) n);  //get the inverse (compute it only once)

  //STEP 1.2: Check valid input size
  if(n != Y.n_elem){  //check if n = Y.n_elem and throw error if not
    throw std::range_error("row number of regressor matrix unequal to number of elements in Y for some partition!");
  }

  //STEP 1.3: Initialize beta using OLS/GRADIENT DESCENT
  if(n>p){ //initialize using OLS or at random if p>=n
    beta = solve(X,Y);  // armadillo function x=solve(A,b) solves the system of linear equations A%*%x = b
    // this is faster than using .i() to get the inverse function of X'X
    //IMPLEMENT C ROUTINE OR GRADIENT DESCENT
    //beta = randn<vec>(p);

  }else{
    beta = randn<vec>(p);   //initialize the coefficients as standard normaXi_Squared_inv, l random variables
    //if the OLS-equation A%*%x = b cannot be solved due to singularity of A=X'X.
  }

  //STEP 1.4 (OPTIONAL): Cut off the beta-vector to zero for values smaller than a certain value
  //                      (This will make the inner loop less expensive, )

  //STEP 1.5: Prepare the objects that will be needed inside of the loop.
  //          Set XXCP, XYCP, inModelBool, inModelIndices, gradients equal to their initial values!
  //          NOTE: -No need to use the 'new' and 'delete' operators here, because this is done within the library
  //                already! So I don't have to worry about it!
  //                -The initialization should be to the correct dimensions in order to tell the compiler which
  //                size of memory it needs to block/use for this object!

  //STEP 1.5.1: define the right dimensions
  mat XXCP(p,p, fill::zeros);
  vec XYCP(p), gradients(p, fill::zeros);     //fill with zeros so we can add multiple times
  Col<unsigned int> inModelIndices(p), helperVec = find(beta); //get all nonzero elements in beta (=1:p unless I cut off the initial OLS estimate)
  Col<unsigned char> inModelBool(p, fill::zeros);


  //STEP 1.5.1: fill them as you should
    //k
  k = helperVec.n_elem;                        //get k at the start of the iterations (will be =n without cutoff)
    //inModelIndices
  inModelIndices.subvec(0,k-1) = helperVec;   //put the first k elements of inModelIndices equal to those non-zero indices
    //inModelBool
  for(unsigned int l=0;l<k;l++){             //set all bools for regressors in the model equal to 1
    inModelBool(inModelIndices(l)) = 1;
  }
    //XYCP
  XYCP = (X.t() * Y) * n_inv; // t(X)%*%Y to get cross-products of X_i and Y, scaled by 1/n
    //XXCP and gradients XXCP part
  for(unsigned int l = 0; l<k ; l++){
    for(unsigned int m = 0; m<k ; m++){
      ind1 = inModelIndices(l); ind2 = inModelIndices(m);       //Get the indices ind1 and ind2 as the 'actual' indices
      if(ind1<=ind2){
        XXCP(ind1,ind2) = dot(X.col(ind1), X.col(ind2))*n_inv;  //only fill the matrix at the positions where
                                                               //we have active regressors at this stage
        gradients(ind1) += XXCP(ind1,ind2);                   //Also: fill the gradient with the X-cross products
      }
    }
  }
    //gradients XYCP and beta part
  gradients = XYCP + beta;    //the parts not involving the sum over non-zero beta crossproducts


  //STEP 2: LOOP AND RESULT COMPUTATION

  //STEP 2.1: loop until convergence criterion achieved
  bool criterion = false;
  unsigned int count = 0;
  while(!criterion && count<=max_iter){
    Rcpp::checkUserInterrupt();                                         //check if User pressed 'Stop'
    beta_update = update_cov(beta, X, Y, lambda, n, p, k,
                             n_inv, XXCP, XYCP, inModelIndices,
                             inModelBool, gradients);                 //get updated beta
    criterion =  ( (max(abs(beta_update - beta))) < precision );       //check if your convergence criterion is satisfied
    count++;                                                        //update the loop counter
    beta = beta_update;  //Update the beta
  }

  //return the beta-vector
  return(beta);
}


//THIS FUNCTION DEPENDA ON ARMADILLO!
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
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_naive(arma::vec &beta, arma::mat &X, arma::vec &Y, double lambda, int n, int p){

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

//THIS FUNCTION DEPENDA ON ARMADILLO!
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
//Coordinate descent on a set
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec coordinateDescent_naive(arma::mat &X, arma::vec &Y, double lambda, double precision, int max_iter){

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

//testing solve from armadillo vs standard R
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec solverTester(arma::mat X, arma::vec Y){
  arma::vec beta = arma::solve(X,Y);
  return(beta);
}

/*** R
timesTwo(42)
*/
