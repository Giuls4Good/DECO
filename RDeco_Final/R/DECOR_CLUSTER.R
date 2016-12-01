#' DECO Clusterized Algorithm (Pure R)
#'
#' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
#' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
#' @param p is the column dimension of X [equivalently, p is the number of regressor variables].
#' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals].
#' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
#' @param r_1 is a tweaking parameter for making the inverse more robust (as we take inverse of XX + r_1*I)
#' @param clust an object obtained by \code{makePSOCKcluster}
#' @param r_2 is a tweaking parameter for making the inverse more robust (as we take inverse of X_MX_M + r_2*I)
#' @param ncores determines the number of threads used on each machine to parallelize computation
#' @param intercept determines whether to include an intercept in the model or not
#' @param refinement determines whether to include the refinement step (Stage 3 of the algorithm)
#' @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
#' @details The algorithm is based on the description in "DECOrrelated feature space partitioning
#'          for distributed sparse regression" in Wang, Dunson, and Leng (2016) if lambda is fixed and
#'          LASSO is used as the penalized regression scheme. The rotated versions of Y and X the authors denote
#'          with Tilde are denoted as X* and Y* in the comments below
#' @note -This implementation uses only R functions.
#'
#' - This implementation is meant to distribute the load of work to several machines. Note that the current implementation
#' does not deal with the problem of storing big matrices; this function is just the starting step and it should be
#' further developed (i.e. reading the matrix chunckwise from a file, C++ implementation,parallelizing on each
#' machine,...).
#' @export
DECO_LASSO_R_CLUSTER <- function(Y, X, p, n, lambda, r_1, clust, r_2 = r_1, ncores=1, intercept=TRUE, refinement = TRUE){
  #***  STEP 1: INITIALIZATION  ***#
  nclusters <- length(clust)

  #**   STEP 1.0 Save given X and Y to compute intercept later on         **#
  X_orig  <- X
  Y_orig  <- Y

  #**   STEP 1.1 Mean Standardization         **#
  Y_stand <- as.vector(Y - mean(Y));
  X <- scale(X,scale=FALSE)[,]

  #**   STEP 1.2 Arbitrary Partitioning       **#
  #X<-X[,sample.int(p, p, replace=FALSE)]  #In case the columns were sorted in some way relevant for
  #predictive power/correlation, randomly disturb the sorting
  #DO NOT USE THIS
  Xi <- mclapply(1:nclusters,
                 function(i){
                   startGroup <- floor(p/m)*(i-1) + 1
                   if(i < m) {
                     endGroup <- floor(p/m)*i
                   } else {
                     endGroup <- p
                   }
                   return( X[,startGroup:endGroup])
                 }, mc.cores=ncores
  )                                                     #Note: matrix(unlist(Partitions), ncol=p) = X

  #**   STEP 1.3 Distribution on cores        **#
  #Only relevant in implementations that actually load the data into different machines

  XiXi_list <- clusterApplyLB(clust, Xi, function(X) {return(X%*%t(X))})

  #**   STEP 2.2 Compute X'X from X(i)'X(i)         **#
  XX <- Reduce( '+', XiXi_list); rm(XiXi_list); gc()
  #**   STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}   **#
  XX_Inverse <- solve(XX + r_1*diag(n)); rm(XX); gc() #obtain the inverse (X'X + rI)^(-1) and save memory
  SVD_Object <- svd(XX_Inverse); rm(XX_Inverse) #save memory  #~PARALLEL~#
  #use SVD and the fact that (X'X + rI) is symmetric (U=V) to get (X'X + rI)^(-0.5).
  #Operation here is equivalent to sqrt(p)*(SVD_Object$u%*%diag(SVD_Object$d)) but faster
  XX_Inverse_Sqrt <- sqrt(p)*t(t(SVD_Object$u)*sqrt(SVD_Object$d))
  XX_Inverse_Sqrt <- XX_Inverse_Sqrt%*%t(SVD_Object$u); rm(SVD_Object);gc()

  #**   STEP 2.4 Compute Y* and X*(i) for each i    **#
  Y <- XX_Inverse_Sqrt%*%Y_stand
  Xi2 <- mclapply(1:nclusters,
                  function(i){
                    return(XX_Inverse_Sqrt%*%(Xi[[i]]))        #Compute XX_Inverse_Sqrt%*%Xi for all i
                  }, mc.cores=ncores
  )
  Xi <- Xi2; rm(Xi2); gc() #store new Xis in old structure and delete the new Xi2 one

  #***  STEP 3: ESTIMATION      ***#

  #**   STEP 3.1 LASSO for coefs on each partition i  **#
  #Prepare data to send
  clust_data <- vector("list", nclusters)
  for (i in 1:nclusters) {
    clust_data[[i]] <- list(X=Xi[[i]], Y=Y, alpha = 1, nlambda = 1, lambda=2*lambda, intercept=FALSE)
  }
  #Compute lasso on the machines
  coefs_list <- clusterApplyLB(clust, clust_data, function(data) {
    require(glmnet)
    myCoefs <- as.vector(coef(glmnet(data$X, data$Y, alpha=data$alpha, nlambda=data$nlambda,
                                 lambda=data$lambda, intercept=data$intercept)))
    return(myCoefs[-1]) #no intercept included
  })
  #**   STEP 3.2 Put together all estimated coefs     **#
  coefs <- unlist(coefs_list)

  #**   STEP 3.3 Insert Intercept into coefs         **#
  if(intercept) {
    coef0 <- mean(Y_orig) - colMeans(X_orig)%*%coefs
    # Insert Intercept
    coefs <- c(coef0, coefs)
  }

  #***  STEP 4: REFINEMENT      ***#
  if (refinement){

    if(intercept){
      #Update the X matrix to inlcude a column for the intercept
      X = cbind(rep(1,n),X)
    }

    #**   STEP 4.1 Check if n<=#nonzero coefs and perform LASSO if so **#
    M = which(abs(coefs) > 0.001 )
    X_M = X[,M]
    if(length(M) >= n){
      coefs[M] = coef(glmnet(X_M, Y_stand, alpha = 1, nlambda = 1, lambda = 2*lambda, intercept = FALSE))[-1]

      M_new = which( abs(coefs[M]) > 0.001 )
      X_M = X_M[,M_new]

      #This can potentially be made faster by doing: M = intersect(M,M_new) or something of this sort! Probably not the bottleneck though
      M = which(abs(coefs) > 0.001 )
    }

    if(length(M) == 0){
      print("All coefficients estimated as 0")
      return(rep(0, length(coefs)))
    }
    #**   STEP 4.2 Run Ridge regression on all non-zero coef vars     **#

    #Apply Ridge Regression to give an updated and hopefully better estimate of the coefficient vector
    coefs[M] = solve(t(X_M) %*% X_M + r_2 * diag(length(M))) %*% t(X_M) %*% Y_stand
  }

  return(coefs)
}
