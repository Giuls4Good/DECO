% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DECO_LASSO_C}
\alias{DECO_LASSO_C}
\title{DECO Parallelized Algorithm (Pure C)}
\usage{
DECO_LASSO_C(Y, X, p, n, lambda, r, ncores = 1L, intercept = TRUE)
}
\description{
This function is deprecated. Use \code{DECO_LASSO_C_PARALLEL} function.
}
\details{
This function is equivalent to \code{DECO_LASSO_C_PARALLEL} function when fixing \code{m=1, ncores=1}.
}

