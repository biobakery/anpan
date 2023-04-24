data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  vector[N] offset;
  matrix[N, N] Lcov;  // cholesky factor of known covariance matrix
  real<lower=0> int_prior_scale;
  vector<lower=0>[K-1] beta_sd;
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
  real centered_cov_intercept;
}
model {
  // likelihood
  target += bernoulli_logit_glm_lpmf(Y | Xc, centered_cov_intercept + offset, beta);

  // priors
  // target += std_normal_lpdf(centered_cov_intercept);
  target += normal_lpdf(centered_cov_intercept | 0, int_prior_scale);

  target += normal_lpdf(beta | 0, beta_sd);
}
generated quantities {
  // actual population-level intercept
  real intercept = centered_cov_intercept - dot_product(means_X, beta);
  vector[N] log_lik;

  vector[N] lin_pred;
  lin_pred = Xc * beta + centered_cov_intercept + offset;

  for (i in 1:N){
    log_lik[i] = bernoulli_logit_glm_lpmf(Y[i] | to_matrix(Xc[i]), centered_cov_intercept + offset[i], beta);
  }
}
