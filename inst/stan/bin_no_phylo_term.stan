data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  matrix[N, N] Lcov;  // cholesky factor of known covariance matrix
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
  vector[Kc] b;  // population-level effects
  real Intercept;
}
model {
  // likelihood
  target += bernoulli_logit_glm_lpmf(Y | Xc, Intercept, b);

  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = bernoulli_logit_glm_lpmf(Y[i] | to_matrix(Xc[i]), Intercept, b);
  }
}
