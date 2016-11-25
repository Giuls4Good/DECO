library(parallel)
clust <- makePSOCKcluster(c("greywagtail",
                             "greyheron",
                             "greypartridge",
                             "greyplover"))
x <- list(X=5,Y=4,Z=8,T=9)
lambda <- clusterApplyLB(clust, x, sqrt)
stopCluster(clust)
