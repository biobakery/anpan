data {
  int<lower=1> N;
  vector[N] pwy_abd;
  real pwy_mean; // observed average log10 pwy abd
  matrix[N, 2] intercept_species; // global intercept + fixed effect of species
  int<lower=1> N_pwy;
  array[N] int<lower=1, upper = N_pwy> pwy_ind;
  array[N] int<lower=0, upper = 1> group_ind; // 0 = ctl, 1 = case
}
transformed data {
  matrix[N, 1] centered_species_abd;
  vector[1] mean_species_abd;

  mean_species_abd[1] = mean(intercept_species[:,2]);
  centered_species_abd[:,1] = intercept_species[:,2] - mean_species_abd[1];
}
parameters {
  vector[1] species_beta;
  real centered_intercept;
  real<lower=0> sigma; // dispersion parameter

  real<lower=0> sd_pwy_int; // group-level standard deviations
  vector[N_pwy] z_pwy_int; // standardized group-level effects
  real<lower=0> sd_pwy_effects; // group-level standard deviations
  vector[N_pwy] z_pwy_effects; // standardized group-level effects
}
transformed parameters {
  vector[N_pwy] pwy_intercepts;
  vector[N_pwy] pwy_effects;
  real lprior = 0;

  pwy_intercepts = sd_pwy_int * z_pwy_int;
  pwy_effects = sd_pwy_effects * z_pwy_effects;

  lprior += student_t_lpdf(centered_intercept | 3, pwy_mean, 2.5);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_pwy_int | 5, 0, 2.5)
            - 1 * student_t_lccdf(0 | 5, 0, 2.5);

  // half std_normal on case effects
  lprior += normal_lpdf(sd_pwy_effects | 0, .333) - 1 * normal_lccdf(0 | 0, .333);
}
model {
  // likelihood
  vector[N] mu = centered_intercept + rep_vector(0.0, N);
  for (i in 1:N) {
    mu[i] += pwy_intercepts[pwy_ind[i]] + group_ind[i] * pwy_effects[pwy_ind[i]];
  }

  target += normal_id_glm_lpdf(pwy_abd | centered_species_abd, mu, species_beta, sigma);

  // priors
  target += lprior;
  target += std_normal_lpdf(z_pwy_int);
  target += std_normal_lpdf(z_pwy_effects);
}
generated quantities {
  // actual population-level intercept
  real global_intercept = centered_intercept - dot_product(mean_species_abd, species_beta);
}
