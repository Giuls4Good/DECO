#ifndef DECO_HEADER
#define DECO_HEADER

arma::vec lassoRCoef(arma::mat&, arma::vec&, double, double, bool, Function&);
arma::vec lassoRCoefParallel(arma::mat*, arma::vec& , double, double, bool, int, Function&);
arma::vec coordinateDescent_naive(arma::mat, arma::vec, double, double, int);
arma::vec refine(arma::vec&, arma::mat&, arma::vec, bool, float, int, int, float);

#endif
