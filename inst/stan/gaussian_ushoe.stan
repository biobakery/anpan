functions {
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
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
parameters {
  vector[K_covariates] b_covariates;  // coefficients for additional (non-horseshoe'd) covariates
  // local parameters for horseshoe prior
  vector[K_genes] zb_genes;
  vector<lower=0>[K_genes] hs_local_genes;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global_genes;  // global shrinkage parameters
  real<lower=0> hs_slab_genes;  // slab regularization parameter
  real<lower=0> sigma;
}
transformed parameters {
  vector[K_genes] b_genes;  // gene coefficients
  b_genes = horseshoe(zb_genes, hs_local_genes, hs_global_genes, hs_scale_slab_genes^2 * hs_slab_genes);
}
model {
  //likelihood
  if (!prior_only) {
    vector[N] covariate_term = X_covariates * b_covariates;
    vector[N] gene_term = X_genes * b_genes;

    vector[N] mu;

    for (n in 1:N) {
      mu[n] = covariate_term[n] + gene_term[n];
    }

    target += normal_lpdf(Y | mu, sigma);
  }

  // priors
  target += normal_lpdf(b_covariates[2:] | 0, 2);
  target += student_t_lpdf(b_covariates[1] | 3, 0, 3.3); //intercept

  target += std_normal_lpdf(zb_genes);
  target += student_t_lpdf(hs_local_genes | hs_df_genes, 0, 1)
    - rows(hs_local_genes) * log(0.5);
  target += student_t_lpdf(hs_global_genes | hs_df_global_genes, 0, hs_scale_global_genes * sigma)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab_genes | 0.5 * hs_df_slab_genes, 0.5 * hs_df_slab_genes);

  target += student_t_lpdf(sigma | 3, 0, 3.3)
    - 1 * student_t_lccdf(0 | 3, 0, 3.3);
}
