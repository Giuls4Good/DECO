lassoCoef <- function(X,Y, nlambda = 1, lambda, intercept=FALSE) {
  myCoefs<- as.vector(coef(glmnet(X, Y, alpha = 1, nlambda = nlambda, lambda=lambda, intercept=intercept)))
  if(intercept) {
    return(myCoefs);
  } else {
    return(myCoefs[-1]);
  }

}

lassoCoefParallel <- function(data, Y, nlambda = 1, lambda, intercept=FALSE, m) {
  coefs <- unlist(mclapply(1:m,
                           function(i){
                             myCoefs <- as.vector(coef(glmnet(data[[i]], Y, alpha = 1, nlambda = nlambda, lambda=lambda, intercept=intercept))) #compute without intercept
                             #2*lambda is used since glmnet is using a slightly different penalizing formula
                             if(intercept) {
                               return(myCoefs)
                             } else {
                               return(myCoefs[-1]) #no intercept included
                             }
                           }, mc.cores=ncores
  ))
  return(coefs)
}
