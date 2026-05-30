#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace Rcpp;
using namespace arma;
using namespace boost::math::quadrature;

// [[Rcpp::export]]
arma::mat arma_exp(arma::mat m) {
  return exp(m);
}

// [[Rcpp::export]]
arma::vec inv_logit(arma::vec x) {
  arma::vec res = 1 / (1 + exp(-x));
  return res;
}

double inv_logit(double x) {
  double res = 1 / (1 + exp(-x));
  return res;
}

double dbern(double yj, double mm) {
  // Simple implementation of dbinom(yj, size = 1, prob = inv_logit(model_mean), log = TRUE)
  return (1-yj) + (-1 + 2*yj)*inv_logit(mm) ;
}

arma::vec dbern(double yj, arma::vec mm) {
  // Simple implementation of dbinom(yj, size = 1, prob = inv_logit(model_mean), log = TRUE)
  return (1-yj) + (-1 + 2*yj)*inv_logit(mm) ;
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

  arma::vec phy_term = -0.5 * pow((phylo_effect_vec - mu_bar_j) / sigma_bar_j, 2.0) -
    log(sigma_bar_j * sqrt2pi);

  arma::vec model_mean = lm_term + phylo_effect_vec;

  arma::vec p = log(dbern(yj, model_mean));
  // This is the binomial llikelihood for a single binary observation, no need to call dbinom.

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
  // Body of this function is the same as above, just a different function
  // signature when I need a single input/output.

  double phy_term = -0.5 * pow((phylo_effect - mu_bar_j) / sigma_bar_j, 2.0) -
    log(sigma_bar_j * sqrt2pi);

  double model_mean = lm_term + phylo_effect;

  double p = log(dbern(yj, model_mean));
  // This is the binomial llikelihood for a single binary observation, no need to call dbinom.

  double res = phy_term + p - offset_term;

  if (1-log_out) {
    res = exp(res);
  }

  return res;
}

// [[Rcpp::export]]
NumericVector bintegrate_finite(double mu_bar_j,
                                double sigma_bar_j,
                                double yj,
                                double lmj,
                                double offset_j,
                                double rt2pi,
                                double lower,
                                double upper) {

  double error;

  auto f1 = [mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi](double pv) {

    // return t*t*1/sqrt(2*M_PI) * std::exp(-t*t / (2));
    return integrand_logistic(pv, mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi, 0);
  };

  double Q = gauss_kronrod<double, 61>::integrate(f1,
                                                  lower,
                                                  upper,
                                                  6,
                                                  1e-9,
                                                  &error);
  // TODO: Decide on:
  // * number of points (trying 61 to start)
  // * max_depth (starting with 6)
  // * tolerance (starting with 1e-8)
  // Error estimates are usually much tighter than the declared limit

  NumericVector res = {Q, error};
  return res;

}

// [[Rcpp::export]]
NumericVector bintegrate_infinite(double mu_bar_j,
                                  double sigma_bar_j,
                                  double yj,
                                  double lmj,
                                  double offset_j,
                                  double rt2pi) {

  double error;

  auto f1 = [mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi](double pv) {

    // return t*t*1/sqrt(2*M_PI) * std::exp(-t*t / (2));
    return integrand_logistic(pv, mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi, 0);
  };

  double Q = gauss_kronrod<double, 61>::integrate(f1,
                                                  -std::numeric_limits<double>::infinity(),
                                                  std::numeric_limits<double>::infinity(),
                                                  6,
                                                  1e-9,
                                                  &error);
  // TODO: Decide on:
  // * number of points (trying 61 to start)
  // * max_depth (starting with 6)
  // * tolerance (starting with 1e-8)
  // Error estimates are usually much tighter than the declared limit

  NumericVector res = {Q, error};
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
                     double rt2pi) {

  // double error;
  // auto f1 = [s](double t) { return t*t*1/sqrt(2*M_PI*s*s) * std::exp(-t*t / (2*s*s)); };
  // double Q = gauss_kronrod<double, 15>::integrate(f1, 0, std::numeric_limits<double>::infinity(), 5, 1.22e-4, &error);

  int n = lm_mean.n_elem;

  double sp2 = pow(sigma_phylo, 2.0);

  arma::vec res(n);
  //^ This will contain the integrated log-likelihood for each observation
  //(indexed by j) with the inputs for the current MCMC iteration i.


  arma::vec mu_bar_v(n);
  // vector of mu_bar_j for each observation

  arma::vec var_v(n);
  // same for var_j

  for (int j=0; j<n; ++j) {

    arma::rowvec sigma_inv_row = sigma_inv_mat.row(j);
    // arma::vec eff_mj = phylo_effects;
    // eff_mj.shed_row(j);
    //There's probably a more efficient way to do this.

    arma::vec mul_res = sigma_inv_row.head_cols(j) * phylo_effects.head(j) +
                        sigma_inv_row.tail_cols(n-j-1) * phylo_effects.tail(n-j-1);

    mu_bar_v(j) = mul_res(0,0);

    arma::vec inv_cor = (sigma_inv_row * cor21_mat.col(j));

    var_v(j) = sp2 - (sp2 * inv_cor(0,0));

  }

  double lmj = 0;
  double yj = 0;
  double mean_j = 0;
  double pv_j = 0;
  double dcov = 0;
  double mu_bar_j = 0;
  double var_j = 0;
  double res_j = 0;
  double sigma_bar_j = 0;
  double error = 0;
  NumericVector int_res(2);
  int re_int = 0;
  arma::vec test_pts(3, fill::ones);
  arma::vec test_vls(3, fill::ones);

  // List int_res;
  // arma::vec mean_vec
  // double offset_j = 0;

  for (int j=0; j<n; ++j) {
    lmj = lm_mean[j];
    yj = Y[j];
    mean_j = effect_means[j];
    pv_j = phylo_effects[j];
    dcov = sp2;
    mu_bar_j = mu_bar_v(j);
    var_j = var_v(j);

    if (var_j <= 0 ) {
      // This happens with off diagonal correlations == 1, or from woodbury numerical fuzz. In that
      // case, the conditional normal likelihood on the phylogenetic effect essentially becomes a dirac
      // delta function, which when integrated yields only the binomial likelihood at mu_bar_j. Inspect
      // vec_integrand_logistic closely if this doesn't make sense.

      res_j = log(dbern(yj, lmj + mu_bar_j));
      res(j) = res_j;
      continue;
    }

    sigma_bar_j = sqrt(var_j);

    // Evaluate an integration offset for stability. Evaluate the integrand at
    // the sampled value of phylo_effect_j on the log scale:
    double offset_j = integrand_logistic(pv_j,
                                         mu_bar_j,
                                         sigma_bar_j,
                                         yj,
                                         lmj,
                                         0,
                                         rt2pi,
                                         1);

    int offset_check = abs(offset_j) > 25;

    if (1-offset_check) {

      int_res = bintegrate_infinite(mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi);

    } else {
      int_res = {1,0};
    }

    re_int = offset_check | (int_res[1] > 1e-3) | int_res[0] < 1e-4 | int_res[0] > 100 ;

    if (re_int) {
      // If there's reason, do some more thoughtful prep work and re-integrate.
      // Maybe the offset was huge, maybe the error was large, etc.

      // Evaluate the integrand at 3 test points around the current value (which
      // is plausible (i.e. not pathological) by virtue of it being an MCMC
      // sample). The likelihood will be fairly close to a parabola, so solve
      // the linear system for those three points. From the parabola coeffcients
      // you can get the argmax and the curvature. Set the integral to be
      // finite, centered on the max, and go out +/- enough for the parabola to
      // go down by 16. exp(-20) ~= 2e-9, so the integral will be very close.
      test_pts[0] = pv_j - .1*sigma_phylo;
      test_pts[1] = pv_j;
      test_pts[2] = pv_j + .1*sigma_phylo;

      test_vls = vec_integrand_logistic(test_pts, mu_bar_j, sigma_bar_j, yj, lmj, 0, rt2pi, 1);

      arma::vec A1 = pow(test_pts, 2.0);
      arma::vec A3(3, fill::ones);
      arma::mat A = join_rows(A1, test_pts, A3);

      arma::vec abc = solve(A,test_vls);

      double a = abc[0];
      double b = abc[1];
      double c = abc[2];

      // https://en.wikipedia.org/wiki/Parabola#As_a_graph_of_a_function :)
      double cent = -b / (2*a);
      double maxv = (4*a*c - b*b) / (4*a);
      offset_j = maxv;
      double d = abs((sqrt(b*b - 4*a*(c - maxv + 20)) - 2*a*cent - b) / (2*a));
      double lo = cent - d;
      double hi = cent + d;

      int_res = bintegrate_finite(mu_bar_j, sigma_bar_j, yj, lmj, offset_j, rt2pi, lo, hi);

    }

    res(j) = log(int_res[0]) + offset_j;

  }

  return res;
}

