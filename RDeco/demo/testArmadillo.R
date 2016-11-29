##TEST ARMADILLO FUNCTIONS
library(expm)
library(rbenchmark)

#STANDARDIZE VECTOR
v <- 1:5000000
stand_v <- as.vector(v-mean(v))
if (all.equal(stand_v,as.vector(standardizeVector(v)))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(Rway <- v-mean(v),standardizeVector(v)) #R version is faster

#STANDARDIZE MATRIX
n<-1000
p <- 5000
M <- matrix(rnorm(n*p,10,5), nrow=n)
stand_m <- scale(M,scale=FALSE)[,]
if (all.equal(stand_m,standardizeMatrix(M))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(scale(M,scale=FALSE),standardizeMatrix(M)) #C version is faster (x2)


#TRANSPOSE MATRIX
M <- matrix(1:((1000)^2), nrow=1000)
if (all.equal(t(M),tMatrix(M))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(t(M),tMatrix(M)) #R version is faster

#MULTIPLY MATRIXCS
A <- matrix(rnorm(1000*500,10,5), nrow=1000, ncol=500)
B <- matrix(rnorm(1000*500,10,5), nrow=500, ncol=1000)
if (all.equal(A%*%B,mulMatrices(A,B))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(A%*%B,mulMatrices(A,B)) #Same speed

#INVERSE MATRIX
M <- matrix(rnorm(1000^2,10,5), nrow=1000)
M_symm <- M%*%t(M)
if (all.equal(solve(M_symm),invSymmMatrix(M_symm))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(solve(M_symm),invSymmMatrix(M_symm)) #C++ version faster (2.5x)

#SQUARE ROOT OF SYMMETRIC MATRIX
n <- 100
A <- matrix(rnorm(n*n,mean=10,sd=5),nrow=n)
A_symm <- A%*%t(A)
if (all.equal(sqrtm(A_symm),squareRootSymmetric(A_symm))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with the square root of a symmetric matrix...")
}
benchmark(sqrtm(A_symm),squareRootSymmetric(A_symm)) #C++ version faster (100x)
