library(rbenchmark)
#test data
n<-5
p<-10
m<-1
sigma<-1
lambda<-0.03
r<-0.01
set.seed(100)
X<-matrix(rnorm(n*p, sd=1), nrow=n)
set.seed(1000)
eps<-rnorm(n, mean = 0, sd = sigma)
set.seed(5000)
activeCoefs<-sample(c(0,1),p, prob=c(0.9, 0.1), replace=TRUE)  ##get coefficients that impact Y
set.seed(10000)
coefs<-activeCoefs*1
Y<-X%*%coefs + eps

res<-DECO_LASSO_R(Y, X, p, n, m, lambda, r, ncores = 8)
res_C<-DECO_LASSO(Y, X, p, n, m, lambda, r, ncores = 8)

#Benchmark
benchmark(DECO_LASSO_R(Y, X, p, n, m, lambda, r, ncores = 8), DECO_LASSO(Y, X, p, n, m, lambda, r, ncores = 8)) #1.3x faster
