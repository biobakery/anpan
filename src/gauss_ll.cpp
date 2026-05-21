#include <Rcpp.h>
#include <RcppArmadillo/>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector llij_gauss_inputs(int p,
                       NumericVector lm_means,
                       NumericMatrix sigma_inv,
                       double sigma_phylo,
                       NumericMatrix cor_arr,
                       NumericVector phylo_effects,
                       double sigma_resid,
                       NumericVector Y,
                       NumericVector cov_mat_diag) {
  // This function aims to replace the block setting up j_df in the original
  // implementation. The lapplys there are the main performance bottleneck of
  // the log likelihood calculations.

  // List res(5);
  NumericVector m1 (p);
  NumericVector s1 (p);
  NumericVector l (p);
  NumericVector s2 (p, sigma_resid);

  // res[1] = m1;
  // res[2] = s1;
  // res[3] = l;
  // res[4] = Y;
  // res[5] = s2;

  int j = 1;

  NumericMatrix::Column inv_col = sigma_inv(_,j);
  NumericMatrix emj = phylo_effects[-j];
  NumericVector res = inv_col * emj;

  // for (int j=0; j<p; ++j) {
  //   NumericMatrix::Column inv_col = sigma_inv(_,j);
  //   NumericMatrix emj = phylo_effects[-j];
  //   m1[j] = (inv_col * emj)(1,1);
  //
  // }
  //
  // List res = List::create(m1, s1, l, Y, s2);

  return res;
}

// [[Rcpp::export]]
NumericVector llij_gauss(NumericVector m1, NumericVector s1, NumericVector l, NumericVector yj, NumericVector s2) {

  int n = m1.length();

  NumericVector div_sq_s1 = 1 / pow(s1, 2);

  NumericVector div_sq_s2 = 1 / pow(s2, 2);

  NumericVector a = -.5 * (div_sq_s1 + div_sq_s2);

  NumericVector b = (m1 * div_sq_s1) + ((yj - l) * div_sq_s2);

  NumericVector c = 0.5 * (-(pow(m1,2) * div_sq_s1) -
      (pow(l,2) - 2*l*yj + pow(yj,2)) * div_sq_s2 -
      (log(2 * M_PI) + 2*log(s1)) -
      (log(2 * M_PI) + 2*log(s2)));

  // NumericVector res = log(sqrt(M_PI) * exp(c - (pow(b, 2) / (4*a))) / sqrt(-a));
  NumericVector res = log(sqrt(M_PI)) + c - (pow(b, 2) / (4*a)) - (log(-a)/2);

  return res;
}

// [[Rcpp::export]]
List get_s12x22inv(NumericMatrix x) {
  int p = x.ncol();
  int nr = x.nrow();

  List res(p);

  for (int i=0; i<p; ++i) {
    NumericVector v = x(_,i);
    v.attr("dim") = Dimension(1,nr);
    res[i] = v;
  }
  return res;
}
