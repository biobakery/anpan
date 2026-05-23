#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec inv_logit(arma::vec x) {
  arma::vec res = 1 / (1 + exp(-x));
  return res;
}

double inv_logit(double x) {
  double res = 1 / (1 + exp(-x));
  return res;
}

// [[Rcpp::export]]
arma::vec vec_integrand_logistic(arma::vec phylo_effect_vec,
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
  // This is the binomial likelihood for a single binary observation, no need to call dbinom.

  arma::vec res = phy_term + p - offset_term;

  if (1-log_out) {
    res = exp(res);
  }

  return res;
}

// [[Rcpp::export]]
double integrand_logistic(double phylo_effect,
                          double mu_bar_j,
                          double sigma_bar_j,
                          double yj,
                          double lm_term,
                          double offset_term,
                          double sqrt2pi,
                          int log_out) {

  double phy_term = -0.5 * pow((phylo_effect - mu_bar_j) / sigma_bar_j, 2.0) - log(sigma_bar_j * sqrt2pi);

  double model_mean = lm_term + phylo_effect;

  double p = log((1 - yj) + (-1 + 2*yj) * inv_logit(model_mean));
  // This is the binomial likelihood for a single binary observation, no need to call dbinom.

  double res = phy_term + p - offset_term;

  if (1-log_out) {
    res = exp(res);
  }

  return res;
}

// [[Rcpp::export]]
arma::vec llij_binom(arma::vec lm_mean,
                     arma::mat sigma_inv_mat, // n x (n-1)
                     arma::mat cor21_mat,     // (n-1) x n
                     double sigma_phylo, // no sigma_resid in Binom version
                     arma::vec phylo_effects,
                     arma::vec Y,
                     arma::vec effect_means, // not used in Gaussian version
                     double s2pi) {

  Function integrt("integrate");

  int n = lm_mean.n_elem;
  double sp2 = pow(sigma_phylo, 2.0);
  arma::vec res(n);

  // arma::mat emj(n-1, n);
  arma::vec mu_bar_v(n);
  arma::vec var_v(n);

  for (int j=0; j<n; ++j) {

    arma::rowvec sigma_inv_row = sigma_inv_mat.row(j);
    // arma::vec eff_mj = phylo_effects;
    // eff_mj.shed_row(j);
    //There's probably a more efficient way to do this.

    arma::vec mul_res = sigma_inv_row.head(j) * phylo_effects.head(j) +
                        sigma_inv_row.head(n-j-1) * phylo_effects.head(n-j-1);

    mu_bar_v(j) = mul_res(0,0);

    arma::vec inv_cor = (sigma_inv_row * cor21_mat.col(j));

    var_v(j) = sp2 - (sp2 * inv_cor(0,0));
  }

  double lmj = 0;
  double yj = 0;
  double mean_j = 0;
  double dcov = 0;
  double mu_bar_j = 0;
  double var_j = 0;
  double res_j = 0;
  double sigma_bar_j = 0;
  List int_res;
  // arma::vec mean_vec
  // double offset_j = 0;

  for (int j=0; j<n; ++j) {
    lmj = lm_mean[j];
    yj = Y[j];
    mean_j = effect_means[j];
    dcov = sp2;
    mu_bar_j = mu_bar_v(j);
    var_j = var_v(j);

    if (var_j <= 0 ) {
      // This happens with off diagonal correlations == 1, or from woodbury numerical fuzz. In that
      // case, the conditional normal likelihood on the phylogenetic effect essentially becomes a dirac
      // delta function, which when integrated yields only the binomial likelihood at mu_bar_j. Inspect
      // vec_integrand_logistic closely if this doesn't make sense.

      // arma::vec model_mean(1,1);
      // model_mean(0,0) = lmj + mu_bar_j;
      // arma::vec ilr = inv_logit(model_mean);
      // double ilm = ilr(0,0);

      res_j = log((1-yj) + (-1 + 2*yj)*inv_logit(lmj + mu_bar_j));
      res(j) = res_j;
      continue;
    }

    sigma_bar_j = sqrt(var_j);

    // Evaluate an integration offset for stability:
    double offset_j = integrand_logistic(effect_means(j),
                                         mu_bar_j,
                                         sigma_bar_j,
                                         yj,
                                         lmj,
                                         0,
                                         s2pi,
                                         1);

    // res_j = int_logis();
    // int_res = integrt(vec_integrand_logistic,
    //                   Named("lower") = R_NegInf,
    //                   Named("upper") = R_PosInf,
    //                   Named("mu_bar_j")         = mu_bar_j,
    //                   Named("sigma_bar_j")      = sigma_bar_j,
    //                   Named("yj")               = yj,
    //                   Named("lm_term")          = lm_mean,
    //                   Named("offset_term")      = offset_j(0,0),
    //                   Named("sqrt2pi")          = s2pi,
    //                   Named("log")              = 0);

    // res(j) = offset_j;

    List int_res = integrt(vec_integrand_logistic)
    // res(j) = res_j;
  }

  return res;
}

// double int_logis()
