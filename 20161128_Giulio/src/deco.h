#ifndef DECO_HEADER
#define DECO_HEADER

arma::vec lassoRCoef(arma::mat&, arma::vec&, double, double, bool, Function&);
arma::vec lassoRCoefParallel(arma::mat*, arma::vec& , double, double, bool, int, Function&);

#endif
