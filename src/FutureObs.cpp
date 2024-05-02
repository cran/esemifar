#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec FARIMAfutureObs(arma::vec& obs, const arma::rowvec& ar_inf, const arma::vec& et, const int meanObs) {

  const int n = obs.size();
  const int m = et.size();

  const int l_ar = ar_inf.size();
  const int l = m - 1;


  arma::vec future_obs = arma::join_cols(obs - meanObs, arma::zeros<arma::vec>(m));
  const int k = m + n;

  const arma::rowvec ar = arma::reverse(ar_inf);

  for (int i = 0; i < m; ++i) {

    future_obs.subvec(n + i, n + i) = ar.subvec(l - i, l_ar - 1) * future_obs.subvec(0, n - 1 + i) + et(i);

  }

  const arma::vec out = future_obs.subvec(n, k - 1) + meanObs;

  return out;

}

