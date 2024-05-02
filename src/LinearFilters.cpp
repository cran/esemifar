#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec ARinftySHORT(const arma::vec& ar, const arma::rowvec& ma, const int max_i) {

  if (max_i == 0) {

    return arma::vec(1, arma::fill::value(-1.0));

  } else {

    const int l_ar = ar.size();
    int l_ma = ma.size();

    arma::rowvec ma2 = ma;
    if (l_ma == 0) {
      ma2 = arma::zeros<arma::rowvec>(1);
      l_ma = 1;
    }

    const arma::rowvec ma_n = arma::reverse<arma::rowvec>(-ma2);

    int times = max_i - l_ar;

    if (times < 0) {
      times = 0;
    }

    arma::vec arInf = arma::join_cols(ar, arma::zeros<arma::vec>(times));

    arma::colvec coefOut = arma::join_cols(arma::zeros<arma::colvec>(l_ma - 1), arma::ones<arma::colvec>(1), arma::zeros<arma::colvec>(max_i));

    const int lc = coefOut.size();

    for (int i = l_ma; i < max_i + l_ma; ++i) {

        coefOut.subvec(i, i) = ma_n * coefOut.subvec(i - l_ma, i - 1) - arInf(i - l_ma);

      }

    return arma::join_cols(arma::vec(1, arma::fill::value(-1.0)), -coefOut.subvec(l_ma, lc - 1));

  }

}


// [[Rcpp::export]]
arma::vec MAinftySHORT(const arma::rowvec& ar, const arma::vec& ma, const int max_i) {

  if (max_i == 0) {

    return arma::ones<arma::vec>(1);

  } else {

    int l_ar = ar.size();
    const int l_ma = ma.size();

    arma::rowvec ar2 = ar;

    if (l_ar == 0) {
      l_ar = 1;
      ar2 = arma::zeros<arma::rowvec>(1);
    }

    const arma::rowvec ar_n = arma::reverse<arma::rowvec>(ar2);

    int times = max_i - l_ma;

    if (times < 0) {
      times = 0;
    }

    const arma::vec ma_s = arma::join_cols(ma, arma::zeros<arma::vec>(times));

    arma::vec coefOut = arma::join_cols(arma::zeros<arma::vec>(l_ar - 1), arma::ones<arma::vec>(1), arma::zeros<arma::vec>(max_i));

    const int lc = coefOut.size();

    for (int i = l_ar; i < max_i + l_ar; ++i) {

      coefOut.subvec(i, i) = ar_n * coefOut.subvec(i - l_ar, i - 1) + ma_s(i - l_ar);

    }

    return arma::join_cols(arma::ones<arma::vec>(1), coefOut.subvec(l_ar, lc - 1));

  }

}

// [[Rcpp::export]]
arma::vec ARinftyLONG(arma::rowvec& ar, const arma::vec& d_coef) {

  arma::vec d2 = arma::reverse(d_coef);

  int l_ar = ar.size();

  ar.subvec(0, l_ar - 1) = -1.0 * ar.subvec(0, l_ar - 1);
  arma::vec c_out = arma::ones<arma::vec>(l_ar);

  for (int i = 1; i < l_ar; ++i) {

    c_out.subvec(i, i) = ar.subvec(0, i) * d2.subvec(l_ar - 1 - i, l_ar - 1);

  }

  c_out.subvec(0, l_ar - 1) = -1.0 * c_out.subvec(0, l_ar - 1);

  return c_out;

}

// [[Rcpp::export]]
arma::rowvec MAinftyLONG(arma::vec& ma, const arma::vec& d_coef) {

  arma::vec d2 = arma::reverse(d_coef);

  int l_ma = ma.size();

  arma::rowvec c_out = arma::ones<arma::rowvec>(l_ma);

  for (int i = 1; i < l_ma; ++i) {

    c_out.subvec(i, i) = ma(i) - (c_out.subvec(0, i - 1) * d2.subvec(l_ma - 1 - i, l_ma - 2));

  }

  return c_out;

}

