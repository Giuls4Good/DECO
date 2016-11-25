#Implementation of DECO using C/C++/R (for a fixed lambda and without refinement step)

#' DECO Parallelized Algorithm
#'
#' This implementation uses C/C++ functions to speed up computation.  Refinement step is missing, lambda is fixed.
#'
#' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
#' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
#' @param p is the column dimension of X [equivalently, p is the number of regressor variables]
#' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals]
#' @param m is the number of groups/blocks you wish to split X into, denoted X(i) for 1 <= i <= m
#' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
#' @param ncores determines the number of cores used on each machine to parallelize computation
#' @param r is a tweaking parameter for making the inverse more robust (as we take inverse of XX + r*I)
#' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
#' @details The algorithm is based on the description in "DECOrrelated feature space partitioning
#'          for distributed sparse regression" in Wang, Dunson, and Leng (2016) if lambda is fixed and
#'          LASSO is used as the penalized regression scheme. The rotated versions of Y and X the authors denote
#'          with Tilde are denoted as X* and Y* in the comments below
#' @export
#' @note -This implementation is suboptimal in that X is already stored in the memory when we start the procedure.
#'       Ideally, one would give in only the LOCATION X is stored at and read it in chunkwise (thus allowing for
#'       larger matrices X, as was intended by the authors).
#'       -The notation #~PARALLEl~# will be introduced in the code whereever one may achieve signficiant gains from
#'       parallelizing
#'       -I could evaluate old expressions in the R version within the mcapply loops! ->saves memory as we write over old
#'       data
#'       -We cannot disturb variable order within the algorithm for output comparison reasons, thus reorder X columns before
#'       running DECO_LASSO (if important)
DECO_LASSO<-function(Y, X, p, n, m, lambda, r, ncores=1, intercept=TRUE){
  #***  STEP 1: INITIALIZATION  ***#

  #**   STEP 1.0 Save given X and Y to compute intercept later on         **#
  if(intercept) {
    X_orig  <- X
    Y_orig  <- Y
  }

  #**   STEP 1.1 Mean Standardization         **#
  Y <- as.vector(Y - mean(Y))
  X <- standardizeMatrix(X)

  #**   STEP 1.2 Arbitrary Partitioning       **#

  #Not really sure if C could bring some improvement here
  Xi<-mclapply(1:m,
               function(i){
                 startGroup<-as.integer(round((p/m)*(i-1))) + 1  #start one position advanced to last endGroup computed
                 endGroup<-as.integer(min(p, round((p/m)*i)))    #make sure you don't have index>vectorsize
                 return( X[,startGroup:endGroup])
               }, mc.cores=ncores
  )                                                     #Note: matrix(unlist(Partitions), ncol=p) = X

  #**   STEP 1.3 Distribution on cores        **#
  #Only relevant in implementations that actually load the data into different machines

  #***  STEP 2: DECORRELATION   ***#

  #**   STEP 2.1 Compute X(i)'X(i) for each i       **#
  XiXi_list<-mclapply(1:m,
                      function(i){
                        return(Xi[[i]]%*%t(Xi[[i]])) #Compute X(i)X(i)' for all i
                      }, mc.cores=ncores
  )

  #**   STEP 2.2 Compute X'X from X(i)'X(i)         **#
  XX<-Reduce( '+', XiXi_list); rm(XiXi_list);gc() #save memory #~PARALLEL~ (NOT REALLY: since sum is a pretty
  #cheap operation, paralellizing and copy all the matrixes can take more time.
  #See: https://stackoverflow.com/questions/39360438/parallel-summation-of-matrices-or-rasters-in-r)

  #**   STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}   **#
  XX_Inverse <- invSymmMatrix(XX + r*diag(n)); rm(XX);gc()
  XX_Inverse_Sqrt <- sqrt(p)*squareRootSymmetric(XX_Inverse)

  #**   STEP 2.4 Compute Y* and X*(i) for each i    **#
  Y<-XX_Inverse_Sqrt%*%Y
  Xi_new<-mclapply(1:m,
                   function(i){
                     return(XX_Inverse_Sqrt%*%(Xi[[i]]))   #Compute XX_Inverse_Sqrt%*%Xi for all i
                   }, mc.cores=ncores
  )
  rm(Xi); gc() #store new Xis in old structure and delete the new Xi2 one

  #***  STEP 3: ESTIMATION      ***#

  #**   STEP 3.1 LASSO for coefs on each partition i  **#
  #**   STEP 3.2 Put together all estimated coefs     **#
  coefs<-unlist(mclapply(1:m,
                         function(i){
                           myCoefs<-as.vector(coef(glmnet(Xi_new[[i]], Y, alpha = 1, nlambda = 1, lambda=2*lambda, intercept=FALSE))) #compute without intercept
                           return(myCoefs[-1]) #no intercept included
                         }, mc.cores=ncores
  ))

  #**   STEP 3.3 Compute Intercept from coefs         **#
  if(intercept) {
    coef0<- mean(Y_orig) - colMeans(X_orig)%*%coefs
    coefs<-c(coef0, coefs)
  }

  #***  STEP 4: REFINEMENT      ***#

  #**   STEP 4.1 Check if n<=#nonzero coefs and perform LASSO if so **#

  #**   STEP 4.2 Run Ridge regression on all non-zero coef vars     **#

  return(coefs)
}
