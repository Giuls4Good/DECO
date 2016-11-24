##TEST ARMADILLO FUNCTIONS
library(expm)
library(rbenchmark)

#STANDARDIZE VECTOR
v <- 1:50000
stand_v <- as.vector(v-mean(v))
if (all.equal(stand_v,as.vector(standardizeVector(v)))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(Rway <- v-mean(v),standardizeVector(v))

#STANDARDIZE MATRIX
M <- matrix(1:10, nrow=10)
stand_m <- M-colMeans(M)
if (all.equal(stand_m,standardizeMatrix(M))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(Rway <- M-colMeans(M),standardizeMatrix(M))


#TRANSPOSE MATRIX
M <- matrix(1:((1000)^2), nrow=1000)
if (all.equal(t(M),tMatrix(M))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(t(M),tMatrix(M))

#MULTIPLY MATRIXCS
A <- matrix(rnorm(1000*5,10,5), nrow=1000, ncol=5)
B <- matrix(rnorm(1000*5,10,5), nrow=5, ncol=1000)
if (all.equal(A%*%B,mulMatrices(A,B))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(A%*%B,mulMatrices(A,B))

#INVERSE MATRIX
M <- matrix(rnorm(100^2,10,5), nrow=100)
M_symm <- M%*%t(M)
if (all.equal(solve(M_symm),invSymmMatrix(M_symm))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with standardize vector...")
}
benchmark(solve(M_symm),invSymmMatrix(M_symm))

#SQUARE ROOT OF SYMMETRIC MATRIX
n <- 50
A <- matrix(rnorm(n*n,mean=10,sd=5),nrow=n)
A_symm <- A%*%t(A)
if (all.equal(sqrtm(A_symm),squareRootSymmetric(A_symm))) {
  print("SUCCESS!")
} else {
  print("ERROR! Something is wrong with the square root of a symmetric matrix...")
}
benchmark(sqrtm(A_symm),squareRootSymmetric(A_symm))
