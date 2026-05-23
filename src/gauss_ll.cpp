#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec llij_gauss_eval(arma::vec m1, arma::vec s1, arma::vec l, arma::vec yj, double s2) {

  arma::vec div_sq_s1 = 1 / pow(s1, 2);

  double div_sq_s2 = 1 / pow(s2, 2);

  arma::vec a = -.5 * (div_sq_s1 + div_sq_s2);

  arma::vec b = (m1 % div_sq_s1) + (div_sq_s2*(yj - l));

  arma::vec c = 0.5 * (-(pow(m1,2) % div_sq_s1) -
                        div_sq_s2*(pow(l,2) - 2*l%yj + pow(yj,2)) -
                        (log(2 * M_PI) + 2*log(s1)) -
                        (log(2 * M_PI) + 2*log(s2)));

  // arma::vec res = log(sqrt(M_PI) * exp(c - (pow(b, 2) / (4*a))) / sqrt(-a));
  arma::vec res = log(sqrt(M_PI)) + c - (pow(b, 2) / (4*a)) - (log(-a)/2);

  return res;
}

// [[Rcpp::export]]
arma::vec llij_gauss(int p,
                   arma::vec lm_means,
                   arma::mat sigma_inv,
                   double sigma_phylo,
                   arma::mat cor_arr,
                   arma::vec phylo_effects,
                   double sigma_resid,
                   arma::vec Y) {
  // This function aims to replace the block setting up j_df in the original
  // implementation. The lapplys there are the main performance bottleneck of
  // the log likelihood calculations.

  arma::vec m1(p);
  arma::vec s1(p);
  // arma::vec s2(p, fill::value(sigma_resid));
  double sp2 = pow(sigma_phylo, 2.0);

  for (int j=0; j<p; ++j) {
    arma::rowvec inv_row = sigma_inv.row(j);

    // vec emj = (phylo_effects); // There's probably a better way to to this with a discontiguous view
    // emj.shed_row(j);
    // vec mean_mul = inv_row * emj;
    arma::vec mean_mul = (inv_row.head(j) * phylo_effects.head(j)) +
      (inv_row.tail(p-j-1) * phylo_effects.tail(p-j-1));

    m1[j] = mean_mul(0,0);

    // vec ivec = join_cols(linspace(0,j-1),linspace(j+1,p-1));
    // vec mean_mul = inv_row * phylo_effects.elem(ivec);
    // vec mean_mul = inv_row * join_cols(phylo_effects.head(j),
    //                                    phylo_effects.tail(p-j-1));

    arma::vec cor_col = cor_arr.col(j);
    arma::vec inv_cor = inv_row * cor_col;

    s1[j] = sqrt(sp2 - (sp2 * inv_cor(0,0)));

  }

  arma::vec res = llij_gauss_eval(m1, s1, lm_means, Y, sigma_resid);

  return res;
}
