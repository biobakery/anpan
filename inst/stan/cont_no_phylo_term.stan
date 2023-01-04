data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  real int_mean;
  real<lower=0> resid_scale;
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
}
model {
  // likelihood
  target += normal_id_glm_lpdf(Y | Xc, centered_cov_intercept, beta, sigma_resid);

  // priors
  target += normal_lpdf(centered_cov_intercept | int_mean, int_prior_scale);
  target += normal_lpdf(beta | 0, beta_sd);

  target += student_t_lpdf(sigma_resid | 5, 0, resid_scale)
    - 1 * student_t_lccdf(0 | 5, 0, resid_scale);
}
generated quantities {
  // actual population-level intercept
  real intercept = centered_cov_intercept - dot_product(means_X, beta);
  vector[N] log_lik;
  vector[N] lin_pred;

  lin_pred = Xc * beta + centered_cov_intercept;

  for (i in 1:N){
    log_lik[i] = normal_id_glm_lpdf(Y[i] | to_matrix(Xc[i]), centered_cov_intercept, beta, sigma_resid);
  }
}
