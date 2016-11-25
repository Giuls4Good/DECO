// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// squareRootSymmetric
arma::mat squareRootSymmetric(arma::mat M);
RcppExport SEXP RDeco_squareRootSymmetric(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(squareRootSymmetric(M));
    return rcpp_result_gen;
END_RCPP
}
// standardizeVector
arma::vec standardizeVector(arma::vec V);
RcppExport SEXP RDeco_standardizeVector(SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(standardizeVector(V));
    return rcpp_result_gen;
END_RCPP
}
// standardizeMatrix
arma::mat standardizeMatrix(arma::mat M);
RcppExport SEXP RDeco_standardizeMatrix(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(standardizeMatrix(M));
    return rcpp_result_gen;
END_RCPP
}
// tMatrix
arma::mat tMatrix(arma::mat M);
RcppExport SEXP RDeco_tMatrix(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(tMatrix(M));
    return rcpp_result_gen;
END_RCPP
}
// mulMatrices
arma::mat mulMatrices(arma::mat A, arma::mat B);
RcppExport SEXP RDeco_mulMatrices(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(mulMatrices(A, B));
    return rcpp_result_gen;
END_RCPP
}
// invSymmMatrix
arma::mat invSymmMatrix(arma::mat M);
RcppExport SEXP RDeco_invSymmMatrix(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(invSymmMatrix(M));
    return rcpp_result_gen;
END_RCPP
}