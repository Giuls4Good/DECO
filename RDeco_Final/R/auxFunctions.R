# Lasso coefficients (singol core implementation)
#
# @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
# @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
# @param nlambda how many lambda should be tested. Right now only \code{nlambda=1} is supported.
# @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
# @param intercept determines whether to include an intercept in the model or not
# @return The coefficients of the Lasso regression.
# @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
# @details This function uses glmnet package to compute Lasso coefficient.
#
# Check \code{lassoCoefParallel} for a parallelized version of this function.
#
# @export
lassoCoef <- function(X,Y, nlambda = 1, lambda, intercept=FALSE) {
  if(nlambda != 1) {
    stop("Currently, only nlambda equal to 1 is supported.")
  }
  myCoefs<- as.vector(coef(glmnet(X, Y, alpha = 1, nlambda = nlambda, lambda=lambda, intercept=intercept)))
  if(intercept) {
    return(myCoefs);
  } else {
    return(myCoefs[-1]);
  }

}

# Lasso coefficients (parallel implementation)
#
# @param data is a list of m matrices X. Each X gives the nxp matrix of regressors,
# each column corresponding to a different regressor
# @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
# @param nlambda how many lambda should be tested. Right now only \code{nlambda=1} is supported.
# @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
# @param intercept determines whether to include an intercept in the model or not
# @return The coefficients of the Lasso regression.
# @author Samuel Davenport, Jack Carter, Giulio Morina, Jeremias Knoblauch
# @details This function uses glmnet package to compute Lasso coefficient in a parallel way.
#
# Check \code{lassoCoef} for a not parallelized version of this function.
#
# @export
lassoCoefParallel <- function(data, Y, nlambda = 1, lambda, intercept=FALSE, m) {
  if(nlambda != 1) {
    stop("Currently, only nlambda equal to 1 is supported.")
  }
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
