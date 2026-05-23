#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec inv_logit(arma::vec x) {
  arma::vec res = 1 / (1 + exp(-x));
  return res;
}

// [[Rcpp::export]]
arma::vec llij_logis_integrand(arma::vec phylo_effect_vec,
                               double mu_bar_j,
                               double sigma_bar_j,
                               double yj,
                               double lm_term,
                               double offset_term,
                               double sqrt2pi,
                               int log_out) {

  arma::vec phy_term = -0.5 * pow((phylo_effect_vec - mu_bar_j) / sigma_bar_j, 2.0) - log(sigma_bar_j * sqrt2pi);

  arma::vec model_mean = lm_term + phylo_effect_vec;

  arma::vec p = log((1 - yj) + (-1 + 2*yj) * inv_logit(model_mean));

  arma::vec res = phy_term + p - offset_term;

  if (1-log_out) {
    res = exp(res);
  }

  return res;
}
