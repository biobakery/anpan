data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // int<lower=1> J_1[N];  // grouping indicator per observation
  matrix[N, N] Lcov;  // cholesky factor of known covariance matrix
  real int_mean;
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
  real<lower=0> sigma;  // dispersion parameter
  real<lower=0> sigma_phylo;  // group-level standard deviations
  vector[N] z_1;
}
transformed parameters {
  vector[N] phylo_effect;  // actual group-level effects
  phylo_effect = (sigma_phylo * (Lcov * z_1));
}
model {
  // likelihood including constants
  vector[N] mu = Intercept + phylo_effect;

  target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);

  // priors including constants
  target += normal_lpdf(Intercept | int_mean, resid_scale);
  target += student_t_lpdf(sigma | 3, 0, resid_scale)
    - 1 * student_t_lccdf(0 | 3, 0, resid_scale);
  target += student_t_lpdf(sigma_phylo | 3, 0, resid_scale)
    - 1 * student_t_lccdf(0 | 3, 0, resid_scale);
  target += std_normal_lpdf(z_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
