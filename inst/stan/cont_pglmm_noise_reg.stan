data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // int<lower=1> J_1[N];  // grouping indicator per observation
  matrix[N, N] Lcov;  // cholesky factor of known covariance matrix
  real int_mean;
  real<lower=0> resid_scale;
  vector[2] reg_gamma_params;
  vector<lower=0>[K-1] beta_sd;
  real<lower=0> int_prior_scale;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] beta;  // population-level effects
  real centered_cov_intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma_resid;  // dispersion parameter
  real<lower=0> sigma_phylo;  // group-level standard deviations
  vector[N] std_phylo_effects;
}
transformed parameters {
  vector[N] phylo_effect;  // actual group-level effects
  phylo_effect = (sigma_phylo * (Lcov * std_phylo_effects));
}
model {
  // likelihood
  vector[N] mu = centered_cov_intercept + phylo_effect;

  target += normal_id_glm_lpdf(Y | Xc, mu, beta, sigma_resid);


  // priors
  target += gamma_lpdf(sigma_phylo / sigma_resid | reg_gamma_params[1], reg_gamma_params[2]);

  target += normal_lpdf(centered_cov_intercept | int_mean, int_prior_scale);
  target += normal_lpdf(beta | 0, beta_sd);

  target += student_t_lpdf(sigma_resid | 5, 0, resid_scale)
    - 1 * student_t_lccdf(0 | 5, 0, resid_scale);

  target += student_t_lpdf(sigma_phylo | 5, 0, resid_scale)
    - 1 * student_t_lccdf(0 | 5, 0, resid_scale);

  target += std_normal_lpdf(std_phylo_effects);
}
generated quantities {
  // actual population-level intercept
  real intercept = centered_cov_intercept - dot_product(means_X, beta);
  array[N] real yrep;
  vector[N] lin_pred;
  for (i in 1:N){
    yrep[i] = normal_rng(centered_cov_intercept + phylo_effect[i] + Xc[i]*beta, sigma_resid);
  }

  lin_pred = Xc * beta + centered_cov_intercept + phylo_effect;
}
