#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec FARIMAforecastAR(arma::vec& obs, const arma::rowvec& ar_inf, const double meanObs, const int m) {

  const int n = obs.size();

  const int l_ar = ar_inf.size();
  const int l = m - 1;

  arma::vec fcasts = arma::join_cols(obs - meanObs, arma::zeros<arma::vec>(m));
  const int k = m + n;

  const arma::rowvec ar = arma::reverse(ar_inf);

  for (int i = 0; i < m; ++i) {

    fcasts.subvec(n + i, n + i) = ar.subvec(l - i, l_ar - 1) * fcasts.subvec(0, n - 1 + i);

  }

  const arma::vec out = fcasts.subvec(n, k - 1) + meanObs;

  return out;

}

