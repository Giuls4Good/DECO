lassoCoef <- function(X,Y, nlambda = 1, lambda, intercept=FALSE) {
  myCoefs<- as.vector(coef(glmnet(X, Y, alpha = 1, nlambda = nlambda, lambda=lambda, intercept=intercept)))
  if(intercept) {
    return(myCoefs);
  } else {
    return(myCoefs[-1]);
  }

}
