functions {
  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slap regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K_covariates;  // number of population-level effects
  matrix[N, K_covariates] X_covariates;  // population-level design matrix
  int<lower=1> K_genes;  // number of population-level effects
  matrix[N, K_genes] X_genes;  // population-level design matrix
  // data for the horseshoe prior
  real<lower=0> hs_df_genes;  // local degrees of freedom
  real<lower=0> hs_df_global_genes;  // global degrees of freedom
  real<lower=0> hs_df_slab_genes;  // slab degrees of freedom
  real<lower=0> hs_scale_global_genes;  // global prior scale
  real<lower=0> hs_scale_slab_genes;  // slab prior scale
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_covariates = K_covariates - 1;
  matrix[N, Kc_covariates] Xc_covariates;  // centered version of X_covariates without an intercept
  vector[Kc_covariates] means_X_covariates;  // column means of X_covariates before centering
  for (i in 2:K_covariates) {
    means_X_covariates[i - 1] = mean(X_covariates[, i]);
    Xc_covariates[, i - 1] = X_covariates[, i] - means_X_covariates[i - 1];
  }
}
parameters {
  vector[Kc_covariates] b_covariates;  // population-level effects
  real Intercept_covariates;  // temporary intercept for centered predictors
  // local parameters for horseshoe prior
  vector[K_genes] zb_genes;
  vector<lower=0>[K_genes] hs_local_genes;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global_genes;  // global shrinkage parameters
  real<lower=0> hs_slab_genes;  // slab regularization parameter
}
transformed parameters {
  vector[K_genes] b_genes;  // population-level effects
  // compute actual regression coefficients
  b_genes = horseshoe(zb_genes, hs_local_genes, hs_global_genes, hs_scale_slab_genes^2 * hs_slab_genes);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_covariates = Intercept_covariates + Xc_covariates * b_covariates;
    // initialize linear predictor term
    vector[N] nlp_genes = X_genes * b_genes;
    // initialize non-linear predictor term
    vector[N] mu;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = nlp_covariates[n] + nlp_genes[n];
    }
    target += bernoulli_logit_lpmf(Y | mu);
  }
  // priors including constants
  target += normal_lpdf(b_covariates | 0, 2);
  target += student_t_lpdf(Intercept_covariates | 3, 0, 2.5);
  target += std_normal_lpdf(zb_genes);
  target += student_t_lpdf(hs_local_genes | hs_df_genes, 0, 1)
    - rows(hs_local_genes) * log(0.5);
  target += student_t_lpdf(hs_global_genes | hs_df_global_genes, 0, hs_scale_global_genes)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab_genes | 0.5 * hs_df_slab_genes, 0.5 * hs_df_slab_genes);
}
generated quantities {
  // actual population-level intercept
  real b_covariates_Intercept = Intercept_covariates - dot_product(means_X_covariates, b_covariates);
}
