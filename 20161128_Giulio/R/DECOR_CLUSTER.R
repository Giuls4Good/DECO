DECO_LASSO_R_CLUSTER <- function(Y, X, p, n, lambda, r_1, clust, r_2 = r_1, ncores=1, intercept=TRUE, refinement = TRUE){
  #***  STEP 1: INITIALIZATION  ***#
  nclusters <- length(clust)

  #**   STEP 1.0 Save given X and Y to compute intercept later on         **#
  X_orig  <- X
  Y_orig  <- Y

  #**   STEP 1.1 Mean Standardization         **#
  Y <- as.vector(Y - mean(Y));
  X <- scale(X,scale=FALSE)[,]

  #**   STEP 1.2 Arbitrary Partitioning       **#
  #X<-X[,sample.int(p, p, replace=FALSE)]  #In case the columns were sorted in some way relevant for
  #predictive power/correlation, randomly disturb the sorting
  #DO NOT USE THIS
  Xi <- mclapply(1:nclusters,
                 function(i){
                   startGroup <- as.integer(round((p/nclusters)*(i-1))) + 1  #start one position advanced to last endGroup computed
                   endGroup <- as.integer(min(p, round((p/nclusters)*i)))    #make sure you don't have index>vectorsize
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
  Y <- XX_Inverse_Sqrt%*%Y
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

    #**   STEP 4.2 Run Ridge regression on all non-zero coef vars     **#
    #Find the indicies of the coefficients that are non-zero
    M = which(coefs != 0);
    #Subset X
    X_M = X[,M]
    #Apply Ridge Regression to give an updated and hopefully better estimate of the coefficient vector
    coefs[M] = solve(t(X_M)%*%X_M + r_2*diag(length(M)))%*% t(X_M) %*% Y
  }

  return(coefs)
}
