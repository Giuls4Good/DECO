pkgname <- "RDeco"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('RDeco')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("invSymmMatrix")
### * invSymmMatrix

flush(stderr()); flush(stdout())

### Name: invSymmMatrix
### Title: Inverse of a matrix
### Aliases: invSymmMatrix

### ** Examples

require(rbenchmark)
M <- matrix(rnorm(1000^2,10,5), nrow=1000)
M_symm <- M%*%t(M)
benchmark(solve(M_symm),invSymmMatrix(M_symm),order='relative')




cleanEx()
nameEx("mulMatrices")
### * mulMatrices

flush(stderr()); flush(stdout())

### Name: mulMatrices
### Title: Multiply two matrices
### Aliases: mulMatrices

### ** Examples

require(rbenchmark)
A <- matrix(rnorm(1000*500,10,5), nrow=1000, ncol=500)
B <- matrix(rnorm(1000*500,10,5), nrow=500, ncol=1000)
benchmark(A%*%B,mulMatrices(A,B),order='relative')




cleanEx()
nameEx("squareRootSymmetric")
### * squareRootSymmetric

flush(stderr()); flush(stdout())

### Name: squareRootSymmetric
### Title: Square root of a symmetric matrix
### Aliases: squareRootSymmetric

### ** Examples

require(expm)
require(rbenchmark)
A <- matrix(rnorm(10000,mean=10,sd=5),nrow=100)
A_symm <- A%*%t(A)
benchmark(sqrtm(A_symm),squareRootSymmetric(A_symm), order='relative')




cleanEx()
nameEx("standardizeMatrix")
### * standardizeMatrix

flush(stderr()); flush(stdout())

### Name: standardizeMatrix
### Title: Standardize a matrix so that its mean is equal to 0
### Aliases: standardizeMatrix

### ** Examples

require(rbenchmark)
M <- matrix(rnorm(1000*5000,10,5), nrow=1000)
benchmark(scale(M,scale=FALSE),standardizeMatrix(M), order='relative')




cleanEx()
nameEx("standardizeVector")
### * standardizeVector

flush(stderr()); flush(stdout())

### Name: standardizeVector
### Title: Standardize a vector
### Aliases: standardizeVector

### ** Examples

require(rbenchmark)
v <- 1:5000000
benchmark({v-mean(v)},standardizeVector(v),scale(v,scale=FALSE), order='relative')




cleanEx()
nameEx("tMatrix")
### * tMatrix

flush(stderr()); flush(stdout())

### Name: tMatrix
### Title: Transpose of a matrix
### Aliases: tMatrix

### ** Examples

require(rbenchmark)
M <- matrix(rnorm(1000*5000,10,5), nrow=1000)
benchmark(t(M),tMatrix(M),order='relative')




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
