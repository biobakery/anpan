data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix

  real<lower=0> resid_scale;
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
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma_resid;  // dispersion parameter
}
model {
  // likelihood including constants
  target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma_resid);

}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = normal_id_glm_lpdf(Y[i] | to_matrix(Xc[i]), Intercept, b, sigma_resid);
  }
}
