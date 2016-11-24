#Implementation 1 of DECO, using R only. Jeremias Knoblauch, 23/11/2016




#' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
#' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
#' @param p is the column dimension of X [equivalently, p is the number of regressor variables]
#' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals]
#' @param m is the number of groups/blocks you wish to split X into, denoted X(i) for 1 <= i <= m
#' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
#' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
#' @details The algorithm is based on the description in "DECOrrelated feature space partitioning 
#'          for distributed sparse regression" in Wang, Dunson, and Leng (2016) if lambda is fixed and
#'          LASSO is used as the penalized regression scheme. The rotated versions of Y and X the authors denote 
#'          with Tilde are denoted as X* and Y* in the comments below
#'          
DECO_LASSO<-function(Y, X, p, n, m, lambda){
  #***  STEP 1: INITIALIZATION  ***#
  
    #**   STEP 1.1 Standardization    **#
  
    #**   STEP 1.2 Partitioning       **#
  
    #**   STEP 1.3 Distribution       **#
  
  #***  STEP 2: DECORRELATION   ***#
    
    #**   STEP 2.1 Compute X(i)'X(i) for each i       **#
  
    #**   STEP 2.2 Compute X'X from X(i)'X(i)         **#
  
    #**   STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}   **#
  
    #**   STEP 2.4 Compute X*(i), Y*(i) for each i    **#
  
  #***  STEP 3: ESTIMATION      ***#
  
    #**   STEP 3.1 LASSO for coefs on each partition i  **#
  
    #**   STEP 3.2 Put together all estimated coefs     **#
  
    #**   STEP 3.3 Compute Intercept from coefs         **#
  
  #***  STEP 4: REFINEMENT      ***#
  
    #**   STEP 4.1 Check if n<=#nonzero coefs and perform LASSO if so **#
  
    #**   STEP 4.2 Run Ridge regression on all non-zero coef vars     **#
}







