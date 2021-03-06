#Implementation 1 of DECO, using R only. Jeremias Knoblauch, 23/11/2016
library(parallel)
library(glmnet)


#' @param Y gives the nx1 vector of observations we wish to approximate with a linear model of type Y = Xb + e
#' @param X gives the nxp matrix of regressors, each column corresponding to a different regressor
#' @param p is the column dimension of X [equivalently, p is the number of regressor variables]
#' @param n is the row dimension of X (and Y) [equivalently, n is the number of observations/individuals]
#' @param m is the number of groups/blocks you wish to split X into, denoted X(i) for 1 <= i <= m
#' @param lambda gives the (fixed) penalty magnitude in the LASSO fit of the algorithm
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
#'       
#'                           
DECO_LASSO<-function(Y, X, p, n, m, lambda, r){
  #***  STEP 1: INITIALIZATION  ***#
    
    #**   STEP 1.1 Mean Standardization         **#
    Y<-as.vector(Y - mean(Y)); X<-(X - colMeans(X))
    
    #**   STEP 1.2 Arbitrary Partitioning       **#
    X<-X[,sample.int(p, p, replace=FALSE)]  #In case the columns were sorted in some way relevant for
                                            #predictive power/correlation, randomly disturb the sorting
    Xi<-mclapply(1:m, 
           function(i){
                startGroup<-as.integer(round((p/m)*(i-1))) + 1  #start one position advanced to last endGroup computed
                endGroup<-as.integer(min(p, round((p/m)*i)))    #make sure you don't have index>vectorsize
                return( X[,startGroup:endGroup])
           }, mc.cores=8
        )                                                     #Note: matrix(unlist(Partitions), ncol=p) = X
    
    #**   STEP 1.3 Distribution on cores        **#
      #Only relevant in implementations that actually load the data into different machines
  
  #***  STEP 2: DECORRELATION   ***#
    
    #**   STEP 2.1 Compute X(i)'X(i) for each i       **#
    XiXi_list<-mclapply(1:m, 
                function(i){
                  return(Xi[[i]]%*%t(Xi[[i]]))        #Compute X(i)X(i)' for all i
                }, mc.cores=8
              ) 

    #**   STEP 2.2 Compute X'X from X(i)'X(i)         **#
    XX<-Reduce( '+', XiXi_list); rm(XiXi_list);gc() #save memory #~PARALLEL~#
    
    #**   STEP 2.3 Use SVD to get (X'X = rI)^{-0.5}   **#
    XX_Inverse<-solve(XX + r*diag(n)); rm(XX);gc() #obtain the inverse (X'X + rI)^(-1) and save memory
    SVD_Object<-svd(XX_Inverse); rm(XX_Inverse) #save memory  #~PARALLEL~#
        #use SVD and the fact that (X'X + rI) is symmetric (U=V) to get (X'X + rI)^(-0.5).
        #Operation here is equivalent to sqrt(p)*(SVD_Object$u%*%diag(SVD_Object$d)) but faster
    XX_Inverse_Sqrt<-sqrt(p)*t(t(SVD_Object$u)*SVD_Object$d); rm(SVD_Object);gc() 
     
    #**   STEP 2.4 Compute Y* and X*(i) for each i    **#
    Y<-(t(XX_Inverse_Sqrt)%*%Y)
    Xi2<-mclapply(1:m, 
                function(i){
                    return(XX_Inverse_Sqrt%*%(Xi[[i]]))        #Compute XX_Inverse_Sqrt%*%Xi for all i
              }, mc.cores=8
        )
    Xi<-Xi2; rm(Xi2); gc() #store new Xis in old structure and delete the new Xi2 one
  
  #***  STEP 3: ESTIMATION      ***#
  
    #**   STEP 3.1 LASSO for coefs on each partition i  **#
    #**   STEP 3.2 Put together all estimated coefs     **#
    coefs<-unlist(mclapply(1:m, 
                    function(i){
                      myCoefs<-as.vector(coef(glmnet(Xi[[i]], Y, alpha = 1, nlambda = 1, lambda=2*lambda, intercept=FALSE))) #compute without intercept
                      return(myCoefs[-1]) #no intercept included
                    }, mc.cores=8
          ))
  
    #**   STEP 3.3 Compute Intercept from coefs         **#
    coef0<- mean(Y) - colMeans(X)%*%coefs
    coefs<-c(coef0, coefs)
  
  #***  STEP 4: REFINEMENT      ***#
  
    #**   STEP 4.1 Check if n<=#nonzero coefs and perform LASSO if so **#
  
    #**   STEP 4.2 Run Ridge regression on all non-zero coef vars     **#
    
  return(coefs)
}


#test data
n<-100
p<-500
m<-1
sigma<-1
lambda<-0.3
r<-0.01
set.seed(100)
X<-matrix(rnorm(n*p, sd=10), nrow=n)
set.seed(1000)
eps<-rnorm(n, mean = 0, sd = sigma)
set.seed(5000)
activeCoefs<-sample(c(0,1),p, prob=c(0.9, 0.1), replace=TRUE)  ##get coefficients that impact Y
set.seed(10000)
coefs<-activeCoefs*rnorm(p)
Y<-X%*%coefs + eps

res<-DECO_LASSO(Y, X, p, n, m, lambda, r)


