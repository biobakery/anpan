#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List llij_gauss_inputs(int p,
                       arma::vec lm_means,
                       arma::mat sigma_inv,
                       double sigma_phylo,
                       arma::mat cor_arr,
                       arma::vec phylo_effects,
                       double sigma_resid,
                       arma::vec Y,
                       arma::vec cov_mat_diag) {
  // This function aims to replace the block setting up j_df in the original
  // implementation. The lapplys there are the main performance bottleneck of
  // the log likelihood calculations.

  arma::vec m1(p);
  arma::vec s1(p);
  arma::vec s2(p, fill::value(sigma_resid));

  for (int j=0; j<p; ++j) {
    arma::rowvec inv_row = sigma_inv.row(j);

    vec emj = (phylo_effects); // There's probably a better way to to this with a discontiguous view
    emj.shed_row(j);
    vec mean_mul = inv_row * emj;
    m1[j] = mean_mul(0,0);

    // vec ivec = join_cols(linspace(0,j-1),linspace(j+1,p-1));
    // vec mean_mul = inv_row * phylo_effects.elem(ivec);
    // vec mean_mul = inv_row * join_cols(phylo_effects.head(j),
    //                                    phylo_effects.tail(p-j-1));

    arma::vec cor_col = cor_arr.col(j);
    arma::vec inv_cor = inv_row * cor_col;

    s1[j] = sqrt(cov_mat_diag(j) - (pow(sigma_phylo,2.0) * inv_cor(0,0)));

  }

  List res = List::create(m1, s1, lm_means, Y, s2);

  return res;
}


// [[Rcpp::export]]
arma::vec llij_gauss_eval(arma::vec m1, arma::vec s1, arma::vec l, arma::vec yj, arma::vec s2) {

  arma::vec div_sq_s1 = 1 / pow(s1, 2);

  arma::vec div_sq_s2 = 1 / pow(s2, 2);

  arma::vec a = -.5 * (div_sq_s1 + div_sq_s2);

  arma::vec b = (m1 % div_sq_s1) + ((yj - l) % div_sq_s2);

  arma::vec c = 0.5 * (-(pow(m1,2) % div_sq_s1) -
                        (pow(l,2) - 2*l%yj + pow(yj,2)) % div_sq_s2 -
                        (log(2 * M_PI) + 2*log(s1)) -
                        (log(2 * M_PI) + 2*log(s2)));

  // arma::vec res = log(sqrt(M_PI) * exp(c - (pow(b, 2) / (4*a))) / sqrt(-a));
  arma::vec res = log(sqrt(M_PI)) + c - (pow(b, 2) / (4*a)) - (log(-a)/2);

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

// [[Rcpp::export]]
arma::vec llij_gauss(int p,
                     arma::vec lm_means,
                     arma::mat sigma_inv,
                     double sigma_phylo,
                     arma::mat cor_arr,
                     arma::vec phylo_effects,
                     double sigma_resid,
                     arma::vec Y,
                     arma::vec cov_mat_diag) {

  List ll_inputs = llij_gauss_inputs(p,
                                     lm_means,
                                     sigma_inv,
                                     sigma_phylo,
                                     cor_arr,
                                     phylo_effects,
                                     sigma_resid,
                                     Y,
                                     cov_mat_diag);

  arma::vec res = llij_gauss_eval(ll_inputs[0], ll_inputs[1], ll_inputs[2], ll_inputs[3], ll_inputs[4]);

  return(res);

}
